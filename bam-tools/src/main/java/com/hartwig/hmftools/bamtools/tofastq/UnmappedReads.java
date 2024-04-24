package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;

import org.apache.logging.log4j.Level;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class UnmappedReads
{
    private static final int NUM_READS_PER_MB = 600;

    private final FastqConfig mConfig;
    private final FastqWriterCache mWriterCache;
    private final FastqWriter mThreadedWriter;

    private final Map<String,SAMRecord> mUnmatchedReads;

    private int mUnmappedReadCount;
    private int mLastCachedCount;

    public UnmappedReads(final FastqConfig config, final FastqWriterCache writerCache)
    {
        mConfig = config;
        mWriterCache = writerCache;

        mUnmatchedReads = Maps.newHashMap();
        mUnmappedReadCount = 0;
        mLastCachedCount = 0;

        mThreadedWriter = mConfig.SplitMode != FileSplitMode.READ_GROUP ? mWriterCache.getThreadedWriter() : null;
    }

    public void run()
    {
        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile)))
        {
            BamSlicer bamSlicer = new BamSlicer(0, true, false, false);
            bamSlicer.setKeepUnmapped();

            // first we find out how many unmapped reads there are, and then calculate how many
            // passes we need
            int numUnmappedReads = 0;
            try(final SAMRecordIterator iterator = samReader.queryUnmapped())
            {
                while(iterator.hasNext())
                {
                    iterator.next();
                    ++numUnmappedReads;
                }
            }
            int numPasses = calcNumPassesRequired(numUnmappedReads);

            // now we do the passes
            for(int i = 0; i < numPasses; ++i)
            {
                try(final SAMRecordIterator iterator = samReader.queryUnmapped())
                {
                    while(iterator.hasNext())
                    {
                        processSamRecord(iterator.next(), numPasses, i);
                    }
                }

                if(!mUnmatchedReads.isEmpty())
                {
                    BT_LOGGER.info("writing {} unmapped & unpaired reads", mUnmatchedReads.size());
                    mUnmatchedReads.values().forEach(mWriterCache::writeUnpairedRead);
                    mUnmatchedReads.clear();
                }

                BT_LOGGER.printf(Level.INFO, "pass #%d, processed %,d unmapped reads", i + 1, mUnmappedReadCount);
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        if(mUnmappedReadCount > 0)
        {
            BT_LOGGER.info("processed {} unmapped reads", mUnmappedReadCount);
        }
    }

    private static final int LOG_COUNT = 1_000_000;

    private void processSamRecord(final SAMRecord read, int numTotalPasses, int passIndex)
    {
        if(read.isSecondaryOrSupplementary())
            return;

        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // unexpected
            return;

        // use the hashcode of the read name to decide if we need to process it now
        // this way we limit the number of reads we process to avoid running out of
        // memory
        if((read.getReadName().hashCode() % numTotalPasses) != passIndex)
        {
            return;
        }

        ++mUnmappedReadCount;
        SAMRecord mate = mUnmatchedReads.remove(read.getReadName());

        if(mate != null)
        {
            mWriterCache.processReadPair(read, mate, mThreadedWriter);
            return;
        }

        mUnmatchedReads.put(read.getReadName(), read);

        int cacheDiff = Math.abs(mLastCachedCount - mUnmatchedReads.size());

        if(cacheDiff >= LOG_COUNT)
        {
            mLastCachedCount = mUnmatchedReads.size();
            BT_LOGGER.printf(Level.DEBUG, "unmapped read cache(%,d), processed(%,d)", mUnmatchedReads.size(), mUnmappedReadCount);
        }
    }

    // in order to avoid running out of memory, we want to limit the number of reads we process each time
    // depending on how much memory we have.
    static int calcNumPassesRequired(int numUnmappedReads)
    {
        int maxMemoryMB = (int)(Runtime.getRuntime().maxMemory() / (1024 * 1024));

        // we need about 1MB per 600 reads
        int numPasses = Math.max(numUnmappedReads / (maxMemoryMB * NUM_READS_PER_MB), 1);

        int numReadsPerPass = numUnmappedReads / numPasses;

        BT_LOGGER.printf(Level.INFO, "found %,d unmapped reads, maxMem(%,dMB), num passes(%d), num reads per pass(%,d)",
                numUnmappedReads, maxMemoryMB, numPasses, numReadsPerPass);

        return numPasses;
    }
}
