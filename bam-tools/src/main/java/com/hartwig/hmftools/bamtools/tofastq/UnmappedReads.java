package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class UnmappedReads
{
    private final FastqConfig mConfig;
    private final FastqWriterCache mWriterCache;
    private final FastqWriter mThreadedWriter;

    private final Map<String,SAMRecord> mUnmatchedReads;

    private int mUnmappedReadCount;

    public UnmappedReads(final FastqConfig config, final FastqWriterCache writerCache)
    {
        mConfig = config;
        mWriterCache = writerCache;

        mUnmatchedReads = Maps.newHashMap();
        mUnmappedReadCount = 0;

        mThreadedWriter = mConfig.SplitMode != FileSplitMode.READ_GROUP ? mWriterCache.getThreadedWriter() : null;
    }

    public void run()
    {
        SamReader samReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        BamSlicer bamSlicer = new BamSlicer(0, true, false, false);
        bamSlicer.setKeepUnmapped();

        try(final SAMRecordIterator iterator = samReader.queryUnmapped())
        {
            while(iterator.hasNext())
            {
                processSamRecord(iterator.next());
            }
        }

        if(!mUnmatchedReads.isEmpty())
        {
            BT_LOGGER.info("writing {} unmapped & unpaired reads", mUnmatchedReads.size());
            mUnmatchedReads.values().forEach(x -> mWriterCache.writeUnpairedRead(x));

            mUnmatchedReads.clear();
        }

        if(mUnmappedReadCount > 0)
        {
            BT_LOGGER.info("processed {} unmapped reads", mUnmappedReadCount);
        }
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag() || read.isSecondaryAlignment())
            return;

        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // unexpected
            return;

        ++mUnmappedReadCount;
        SAMRecord mate = mUnmatchedReads.remove(read.getReadName());

        if(mate != null)
        {
            mWriterCache.processReadPair(read, mate, mThreadedWriter);
            return;
        }

        mUnmatchedReads.put(read.getReadName(), read);
    }
}
