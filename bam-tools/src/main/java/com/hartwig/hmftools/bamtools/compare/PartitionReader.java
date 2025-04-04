package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.CompareUtils.compareReads;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.ORIG_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.VALUE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

public class PartitionReader implements Runnable
{
    private final String mName;
    private final CompareConfig mConfig;
    private final BamPartition mBamPartition;
    final BamReaderProvider mOrigBamReaderProvider;
    final BamReaderProvider mNewBamReaderProvider;

    private final Map<ReadKey, SAMRecord> mOrigBamReads = new HashMap<>();
    private final Map<ReadKey, SAMRecord> mNewBamReads = new HashMap<>();
    private final ReadWriter mReadWriter;

    @Nullable
    private final UnmatchedReadHandler mUnmatchedReadHandler;
    private final boolean mLogReadIds;

    private long mReadCount;
    private long mNextLogReadCount;

    private final Statistics mStats;

    public PartitionReader(
            final String name,
            final CompareConfig config, final BamPartition bamPartition,
            final BamReaderProvider bamAReaderProvider, final BamReaderProvider bamBReaderProvider,
            final ReadWriter readWriter, @Nullable final UnmatchedReadHandler unmatchedReadHandler)
    {
        mName = name;
        mConfig = config;
        mBamPartition = bamPartition;
        mOrigBamReaderProvider = bamAReaderProvider;
        mNewBamReaderProvider = bamBReaderProvider;
        mReadWriter = readWriter;
        mUnmatchedReadHandler = unmatchedReadHandler;
        mStats = new Statistics();

        mReadCount = 0;
        mNextLogReadCount = LOG_COUNT;

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public Statistics stats() { return mStats; }

    public void run()
    {
        // BT_LOGGER.debug("processing {}", mName);

        // we process the records partition by partition
        // reads are stored inside a hash table and looked up by the read id
        try(final SAMRecordIterator origBamItr = mBamPartition.iterator(mOrigBamReaderProvider.getBamReader());
            final SAMRecordIterator newBamItr = mBamPartition.iterator(mNewBamReaderProvider.getBamReader()))
        {
            long origBamReadCount = 0;
            long newBamReadCount = 0;
            int origReadAlignStart = -1;
            int newReadAlignStart = -1;
            while(origBamItr.hasNext() || newBamItr.hasNext())
            {
                // decide which one to advance, or both
                boolean advanceOrig;
                boolean advanceNew;

                if(origReadAlignStart == newReadAlignStart && origBamItr.hasNext() && newBamItr.hasNext())
                {
                    // if the align start are the same try to advance both. This logic is required such that processing
                    // of unmapped region is not massively slowed down by only advancing one side all the time
                    advanceOrig = true;
                    advanceNew = true;
                }
                else
                {
                    // the one with smaller alignment start is advanced first
                    // this helps to process the same genomic location together
                    advanceOrig = !newBamItr.hasNext() || (origBamItr.hasNext() && (origReadAlignStart < newReadAlignStart));
                    advanceNew = !advanceOrig;
                }

                if(advanceOrig)
                {
                    origBamReadCount++;
                    final SAMRecord origBamRead = origBamItr.next();
                    processOrigBamRecord(origBamRead);
                    origReadAlignStart = origBamRead.getAlignmentStart();
                }

                if(advanceNew)
                {
                    newBamReadCount++;
                    final SAMRecord newBamRead = newBamItr.next();
                    processNewBamRecord(newBamRead);
                    newReadAlignStart = newBamRead.getAlignmentStart();
                }

                // check if we need to dump the reads to unmatched read handler, this is to protect against running out of memory
                if(mUnmatchedReadHandler != null && mOrigBamReads.size() + mNewBamReads.size() > mConfig.MaxCachedReadsPerThread)
                {
                    BT_LOGGER.trace("max cache read limit exceeded, writing reads to hash bam");
                    mUnmatchedReadHandler.handleOrigBamReads(mOrigBamReads.values());
                    mUnmatchedReadHandler.handleNewBamReads(mNewBamReads.values());
                    mOrigBamReads.clear();
                    mNewBamReads.clear();
                }
            }
            completePartition(origBamReadCount, newBamReadCount);
        }

        if(mStats.OrigReadCount > 0 || mStats.NewReadCount > 0)
        {
            BT_LOGGER.trace("{} complete: origReads({}) newReads({}) diff({})",
                    mName, mStats.OrigReadCount, mStats.NewReadCount, mStats.DiffCount);
        }
    }

    private void processOrigBamRecord(final SAMRecord origBamRead)
    {
        if(mLogReadIds && mConfig.LogReadIds.contains(origBamRead.getReadName()))
        {
            BT_LOGGER.debug("orig bam: specific readId({})", origBamRead.getReadName());
        }

        if(excludeRead(origBamRead))
            return;

        checkLogReadCounts();

        // create a key
        ReadKey readKey = ReadKey.from(origBamRead);

        // look it up see if can find corresponding beta record
        SAMRecord newBamRead = mNewBamReads.remove(readKey);

        if(newBamRead != null)
        {
            // found matching
            checkReadDetails(origBamRead, newBamRead);
        }
        else
        {
            // add to read map
            mOrigBamReads.put(readKey, origBamRead);
        }
    }

    private void processNewBamRecord(final SAMRecord newBamRead)
    {
        if(mLogReadIds && mConfig.LogReadIds.contains(newBamRead.getReadName()))
        {
            BT_LOGGER.debug("new bam: specific readId({})", newBamRead.getReadName());
        }

        if(excludeRead(newBamRead))
            return;

        checkLogReadCounts();

        // create a key
        ReadKey readKey = ReadKey.from(newBamRead);

        // look it up see if can find corresponding beta record
        SAMRecord origRead = mOrigBamReads.remove(readKey);

        if(origRead != null)
        {
            // found matching
            checkReadDetails(origRead, newBamRead);
        }
        else
        {
            // add to beta read map
            mNewBamReads.put(readKey, newBamRead);
        }
    }

    private void completePartition(long origBamReads, long newBamReads)
    {
        /*BT_LOGGER.printf(Level.DEBUG, "partition(%s) complete, orig bam reads(%,d), new bam reads(%,d)",
                mName, origBamReads, newBamReads);*/

        if(mUnmatchedReadHandler != null)
        {
            // any reads not matched up is palmed off to the unmatched reads handler
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
        if(read.isSecondaryAlignment())
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
                    mBamPartition, mReadCount, mStats.OrigReadCount, mStats.NewReadCount);
        }
    }
}
