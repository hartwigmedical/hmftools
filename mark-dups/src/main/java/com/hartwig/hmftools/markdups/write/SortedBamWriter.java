package com.hartwig.hmftools.markdups.write;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

public class SortedBamWriter
{
    private final SortedBamConfig mConfig;

    private int mCapacity; // only purpose of capacity is to define when the cache is checked
    private int mCapacityGrowThreshold;
    private int mCapacityShrinkThreshold;
    private int mCapacityThresholdCheck;
    private final int mPositionBuffer;

    private final List<SAMRecord> mRecords;
    private final List<SAMRecord> mSortedRecords;
    private final SAMFileWriter mWriter;

    // state to control which reads can be processed
    private String mCurrentChromosome; // reads on another chromosome cannot be
    private int mLastWrittenPosition; // last write from the sorted cache - no lower positions can be processed
    private int mUpperBoundPosition; // most recent read, forming the upper bound for accepting reads
    private int mUpperSortablePosition; // upper position on reads ready for sorting and writing
    private int mMinCachedPosition;

    // stats
    private int mWriteCount;
    private long mReadCount;
    private long mReadsWritten;
    private int mMaxCached;
    private int mMaxWrite;
    private double mAvgCacheSize;
    private long mCacheAssessCount;

    private boolean mPerfDebug;

    public SortedBamWriter(final SortedBamConfig config, @Nullable final SAMFileWriter writer)
    {
        mConfig = config;
        setCapacity(mConfig.Capacity);
        mPositionBuffer = mConfig.PositionBuffer;

        mRecords = Lists.newArrayListWithCapacity(mCapacity);
        mSortedRecords = Lists.newArrayListWithCapacity(mCapacity);
        mWriter = writer;
        mCurrentChromosome = "";
        mLastWrittenPosition = 0;
        mUpperSortablePosition = 0;
        mMinCachedPosition = 0;
        mReadsWritten = 0;
        mWriteCount = 0;
        mReadCount = 0;
        mMaxCached = 0;
        mMaxWrite = 0;
        mAvgCacheSize = 0;
        mCacheAssessCount = 0;
        mPerfDebug = false;
    }

    public void initialiseStartPosition(final String chromosome, int startPosition)
    {
        mCurrentChromosome = chromosome;
        mLastWrittenPosition = 0;
        mMinCachedPosition = startPosition;
        mUpperBoundPosition = startPosition + mPositionBuffer;
        mUpperSortablePosition = startPosition;
    }

    public void setUpperBoundPosition(int position) { mUpperBoundPosition = position; }
    public void setUpperSortablePosition(int position) { mUpperSortablePosition = position - mConfig.ReadPosCacheBuffer; }

    public void togglePerfDebug() { mPerfDebug = true; }

    public boolean canWriteRecord(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() && read.getMateUnmappedFlag())
            return false;

        if(!mCurrentChromosome.equals(read.getReferenceName()))
            return false;

        if(!positionWithin(read.getAlignmentStart(), mLastWrittenPosition, mUpperBoundPosition))
            return false;

        return true;
    }

    private static final int LOG_COUNT = 1000000;

    public void addRecord(final SAMRecord read)
    {
        mRecords.add(read);
        mMaxCached = max(mMaxCached, mRecords.size());

        ++mReadCount;

        if(mPerfDebug && (mReadCount % LOG_COUNT) == 0)
        {
            MD_LOGGER.debug("sorted read cache chr({}:{}) records({}) avgWrite({}) assess(avg={} count={})",
                    mCurrentChromosome, mLastWrittenPosition, mRecords.size(), avgWriteCount(), avgAssessSize(), mCacheAssessCount);
        }

        if(mMinCachedPosition == 0)
            mMinCachedPosition = read.getAlignmentStart();
        else
            mMinCachedPosition = min(mMinCachedPosition, read.getAlignmentStart());

        if(mRecords.size() < mCapacityThresholdCheck)
            return;

        if(mUpperSortablePosition <= mMinCachedPosition)
            return;

        // quick scan to see how many records could be written
        int sortableCount = 0;
        for(SAMRecord cachedRead : mRecords)
        {
            if(cachedRead.getAlignmentStart() <= mUpperSortablePosition)
                ++sortableCount;
        }

        double newCacheAssessCount = mCacheAssessCount + 1;
        double adjustAvgCacheSize = mCacheAssessCount / newCacheAssessCount * mAvgCacheSize;
        mAvgCacheSize = adjustAvgCacheSize + (mRecords.size() / newCacheAssessCount);
        ++mCacheAssessCount;

        if(sortableCount < mConfig.MinSortWriteCount)
            return;

        // gather all reads prior to this position
        int index = 0;
        mMinCachedPosition = mUpperSortablePosition + 1;
        while(index < mRecords.size())
        {
            SAMRecord cachedRead = mRecords.get(index);

            if(cachedRead.getAlignmentStart() <= mUpperSortablePosition)
            {
                mSortedRecords.add(cachedRead);
                mRecords.remove(index);
            }
            else
            {
                // cannot break since the latter records may have earlier positions
                mMinCachedPosition = min(mMinCachedPosition, cachedRead.getAlignmentStart());
                ++index;
            }
        }

        // sort and write them
        sortAndWrite(mSortedRecords);

        checkCapacity();
    }

    private void checkCapacity()
    {
        // check whether capacity needs to grow or shrink
        if(mRecords.size() > mCapacityGrowThreshold)
        {
            MD_LOGGER.debug("sorted read cache chr({}) pos(lastWrite={} upperSort={} upperBound={}) cached({}) growing capacity({} -> {})",
                    mCurrentChromosome, mLastWrittenPosition, mUpperSortablePosition, mUpperBoundPosition,
                    mRecords.size(), mCapacity, mCapacity + mConfig.Capacity);

            setCapacity(mCapacity + mConfig.Capacity);
        }
        else if(mCapacity > mConfig.Capacity && mRecords.size() < mCapacityShrinkThreshold)
        {
            MD_LOGGER.debug("sorted read cache chr({}) pos(lastWrite={} upperSort={} upperBound={}) cached({}) reducing capacity({} -> {})",
                    mCurrentChromosome, mLastWrittenPosition, mUpperSortablePosition, mUpperBoundPosition,
                    mRecords.size(), mCapacity, mConfig.Capacity);

            setCapacity(mConfig.Capacity);
        }
    }

    private void setCapacity(int capacity)
    {
        mCapacity = capacity;
        mCapacityGrowThreshold = (int)(capacity * mConfig.CapacityGrowPercent);
        mCapacityShrinkThreshold = (int)(max(capacity * mConfig.CapacityShrinkPercent, 1));
        mCapacityThresholdCheck = (int)max(capacity * mConfig.CapacityCheckPercent, 1);
    }

    private void sortAndWrite(final List<SAMRecord> records)
    {
        if(records.isEmpty())
            return;

        mReadsWritten += records.size();
        Collections.sort(records, new SAMRecordCoordinateComparator());

        if(records.get(0).getAlignmentStart() < mLastWrittenPosition)
        {
            MD_LOGGER.error("sorted BAM cache({}) writing earlier read({}) from {} records",
                    toString(), readToString(records.get(0)), mRecords.size());
            System.exit(1);
        }

        if(mWriter != null)
            records.forEach(x -> mWriter.addAlignment(x));

        mLastWrittenPosition = records.get(records.size() - 1).getAlignmentStart();
        mMaxWrite = max(mMaxWrite, records.size());

        records.clear();
        ++mWriteCount;
    }

    public void flush()
    {
        // sort and write all cached records
        sortAndWrite(mRecords);
    }

    public long written() { return mReadsWritten; }
    public long writeCount() { return mWriteCount; }
    public int cached() { return mRecords.size(); }
    public int capacity() { return mCapacity; }
    public int maxCache() { return mMaxCached; }
    public int maxWrite() { return mMaxWrite; }
    public int avgAssessSize() { return (int)mAvgCacheSize; }

    public int avgWriteCount() { return mWriteCount > 0 ? (int)(mReadsWritten / (double)mWriteCount) : 0; }

    public String toString()
    {
        return format("chr(%s) lastWritePos(%d) upper(sort=%d bound=%d) cached(%d/%d max=%d) avgWriteCount(%d from %d) thresholds(check=%d grow=%d shrink=%d)",
                mCurrentChromosome, mLastWrittenPosition, mUpperSortablePosition, mUpperBoundPosition, mRecords.size(), mCapacity, mMaxCached,
                avgWriteCount(), mWriteCount, mCapacityThresholdCheck, mCapacityGrowThreshold, mCapacityShrinkThreshold);
    }
}
