package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.CompareUtils.compareReads;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.ORIG_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.VALUE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionReader implements Callable<Void>
{
    private final String mName;
    private final CompareConfig mConfig;

    private final File mOrigBamFile;
    private final File mNewBamFile;
    private final ChrBaseRegion mPartition;

    private final Map<ReadKey,SAMRecord> mOrigBamReads;
    private final Map<ReadKey,SAMRecord> mNewBamReads;
    private final ReadWriter mReadWriter;

    @Nullable
    private final UnmatchedReadHandler mUnmatchedReadHandler;
    private final boolean mLogReadIds;

    private long mReadCount;
    private long mNextLogReadCount;

    private final Statistics mStats;

    protected static final String NAME_UNMAPPED = "fully_unmapped";

    public PartitionReader(
            final String name, final CompareConfig config, @Nullable final ChrBaseRegion bamPartition,
            final File origBamFile, final File newBamFile,
            final ReadWriter readWriter, @Nullable final UnmatchedReadHandler unmatchedReadHandler)
    {
        mName = name;
        mConfig = config;
        mPartition = bamPartition;

        mOrigBamFile = origBamFile;
        mNewBamFile = newBamFile;
        mReadWriter = readWriter;
        mUnmatchedReadHandler = unmatchedReadHandler;
        mStats = new Statistics();

        mOrigBamReads = new HashMap<>();
        mNewBamReads = new HashMap<>();

        mReadCount = 0;
        mNextLogReadCount = LOG_COUNT;

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public Statistics stats() { return mStats; }

    @Override
    public Void call()
    {
        BT_LOGGER.trace("processing {}", mName);

        SamReader origSamReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(mOrigBamFile);

        SamReader newSamReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(mNewBamFile);

        // reads are stored inside a hash table and looked up by the read ID
        SAMRecordIterator origBamIter;
        SAMRecordIterator newBamIter;

        boolean fullyUnmappedReads = mName.equals(NAME_UNMAPPED);

        if(mPartition != null)
        {
            origBamIter = origSamReader.query(mPartition.chromosome(), mPartition.start(), mPartition.end(), false);
            newBamIter = newSamReader.query(mPartition.chromosome(), mPartition.start(), mPartition.end(), false);
        }
        else
        {
            if(fullyUnmappedReads)
            {
                origBamIter = origSamReader.queryUnmapped();
                newBamIter = newSamReader.queryUnmapped();
            }
            else
            {
                // all records (in the hash BAM)
                origBamIter = origSamReader.iterator();
                newBamIter = newSamReader.iterator();
            }
        }

        int origReadAlignStart = -1;
        int newReadAlignStart = -1;

        while(origBamIter.hasNext() || newBamIter.hasNext())
        {
            // decide which one to advance, or both
            boolean advanceOrig;
            boolean advanceNew;

            if(fullyUnmappedReads)
            {
                // while both BAMs have records, move ahead in one by a fixed amount to try to match off as many reads as possible
                // even though they aren't ordered
                if(origBamIter.hasNext() && newBamIter.hasNext())
                {
                    advanceOrig = true;
                    advanceNew = true;
                }
                else
                {
                    advanceOrig = origBamIter.hasNext();
                    advanceNew = !advanceOrig;
                }
            }
            else
            {
                if(origReadAlignStart == newReadAlignStart && origBamIter.hasNext() && newBamIter.hasNext())
                {
                    // if the align start are the same try to advance both. This logic is required such that processing
                    // of unmapped region is not massively slowed down by only advancing one side all the time
                    advanceOrig = true;
                    advanceNew = true;
                }
                else
                {
                    // the one with smaller alignment start is advanced first, so as to process the same genomic location together
                    advanceOrig = !newBamIter.hasNext() || (origBamIter.hasNext() && (origReadAlignStart < newReadAlignStart));
                    advanceNew = !advanceOrig;
                }
            }

            if(advanceOrig)
            {
                SAMRecord origBamRead = origBamIter.next();
                processOrigBamRecord(origBamRead);
                origReadAlignStart = origBamRead.getAlignmentStart();
            }

            if(advanceNew)
            {
                SAMRecord newBamRead = newBamIter.next();
                processNewBamRecord(newBamRead);
                newReadAlignStart = newBamRead.getAlignmentStart();
            }

            // check if we need to dump the reads to unmatched read handler, this is to protect against running out of memory
            int readCacheSize = mOrigBamReads.size() + mNewBamReads.size();

            if(mUnmatchedReadHandler != null && readCacheSize > mConfig.MaxCachedReadsPerThread)
            {
                BT_LOGGER.debug("{}: cached reads({}) exceeds limit, writing reads to hash bam", mName, readCacheSize);
                mUnmatchedReadHandler.handleOrigBamReads(mOrigBamReads.values());
                mUnmatchedReadHandler.handleNewBamReads(mNewBamReads.values());
                mOrigBamReads.clear();
                mNewBamReads.clear();
            }
        }

        completePartition();

        if(mStats.OrigReadCount > 0 || mStats.NewReadCount > 0)
        {
            BT_LOGGER.trace("{} complete: origReads({}) newReads({}) diff({})",
                    mName, mStats.OrigReadCount, mStats.NewReadCount, mStats.DiffCount);
        }

        try
        {
            origSamReader.close();
            newSamReader.close();
        }
        catch(Exception e)
        {
            BT_LOGGER.error("failed to close BAMs");
        }

        return null;
    }

    private void processOrigBamRecord(final SAMRecord origRead)
    {
        if(mPartition != null && !mPartition.containsPosition(origRead.getAlignmentStart()))
            return;

        if(mLogReadIds && mConfig.LogReadIds.contains(origRead.getReadName()))
        {
            BT_LOGGER.debug("orig bam: specific readId({})", origRead.getReadName());
        }

        if(excludeRead(origRead))
            return;

        checkLogReadCounts();

        // create a key
        ReadKey readKey = ReadKey.from(origRead);

        // look it up see if can find corresponding beta record
        SAMRecord newBamRead = mNewBamReads.remove(readKey);

        if(newBamRead != null)
        {
            // found matching
            checkReadDetails(origRead, newBamRead);
        }
        else
        {
            // add to read map
            mOrigBamReads.put(readKey, origRead);
        }
    }

    private void processNewBamRecord(final SAMRecord newRead)
    {
        if(mPartition != null && !mPartition.containsPosition(newRead.getAlignmentStart()))
            return;

        if(mLogReadIds && mConfig.LogReadIds.contains(newRead.getReadName()))
        {
            BT_LOGGER.debug("new bam: specific readId({})", newRead.getReadName());
        }

        if(excludeRead(newRead))
            return;

        checkLogReadCounts();

        // create a key
        ReadKey readKey = ReadKey.from(newRead);

        // look it up see if can find corresponding beta record
        SAMRecord origRead = mOrigBamReads.remove(readKey);

        if(origRead != null)
        {
            // found matching
            checkReadDetails(origRead, newRead);
        }
        else
        {
            // add to beta read map
            mNewBamReads.put(readKey, newRead);
        }
    }

    private void completePartition()
    {
        /*BT_LOGGER.printf(Level.DEBUG, "partition(%s) complete, orig bam reads(%,d), new bam reads(%,d)",
                mName, origBamReads, newBamReads);*/

        if(mUnmatchedReadHandler != null)
        {
            // unmatched reads are then handled by the unmatched reads handler
            if(!mOrigBamReads.isEmpty())
            {
                mUnmatchedReadHandler.handleOrigBamReads(mOrigBamReads.values());
            }
            if(!mNewBamReads.isEmpty())
            {
                mUnmatchedReadHandler.handleNewBamReads(mNewBamReads.values());
            }
        }
        else
        {
            // write them out as mismatches
            mOrigBamReads.values().forEach(read -> mReadWriter.writeComparison(read, ORIG_ONLY, null));
            mStats.OrigReadCount += mOrigBamReads.size();
            mStats.DiffCount += mOrigBamReads.size();

            mNewBamReads.values().forEach(read -> mReadWriter.writeComparison(read, NEW_ONLY, null));
            mStats.NewReadCount += mNewBamReads.size();
            mStats.DiffCount += mNewBamReads.size();
        }

        mOrigBamReads.clear();
        mNewBamReads.clear();
    }

    private void checkReadDetails(final SAMRecord origRead, final SAMRecord newRead)
    {
        ++mStats.OrigReadCount;
        ++mStats.NewReadCount;

        // note that the ReadKey of both reads match
        if(mConfig.IgnoreReduxAlterations && (origRead.hasAttribute(UNMAP_ATTRIBUTE) || newRead.hasAttribute(UNMAP_ATTRIBUTE)))
            return;

        List<String> diffs = compareReads(origRead, newRead, mConfig);

        if(!diffs.isEmpty())
        {
            ++mStats.DiffCount;
            mReadWriter.writeComparison(origRead, VALUE, diffs);
        }
    }

    private boolean excludeRead(final SAMRecord read)
    {
        if(mConfig.IgnoreSecondaryReads && read.isSecondaryAlignment())
            return true;

        if(mConfig.IgnoreReduxAlterations)
        {
            if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE) || read.hasAttribute(UNMAP_ATTRIBUTE))
                return true;
        }

        if(mConfig.IgnoreConsensusReads && read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
            return true;

        if(mConfig.IgnoreSupplementaryReads && read.getSupplementaryAlignmentFlag())
            return true;

        return false;
    }

    private static final int LOG_COUNT = 1_000_000;

    private void checkLogReadCounts()
    {
        ++mReadCount;

        if(mReadCount >= mNextLogReadCount)
        {
            mNextLogReadCount += LOG_COUNT;

            BT_LOGGER.debug("partition({}) reads processed total({}) orig({}) new({})",
                    mPartition, mReadCount, mStats.OrigReadCount, mStats.NewReadCount);
        }
    }
}
