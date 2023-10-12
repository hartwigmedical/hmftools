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
    private final int mCapacityStart;
    private int mCapacity;
    private int mCapacityGrowThreshold;
    private int mCapacityShrinkThreshold;
    private int mCapacityThresholdCheck;
    private final int mPositionBuffer;

    private final List<SAMRecord> mRecords;
    private final List<SAMRecord> mSortedRecords;
    private final SAMFileWriter mWriter;
    private String mCurrentChromosome;
    private int mLastWrittenPosition;
    private int mMinCachedPosition;
    private int mUpperBoundPosition;

    // stats
    private int mWriteCount;
    private long mReadsWritten;
    private int mMaxCached;

    public SortedBamWriter(final int capacity, final int positionBuffer, @Nullable final SAMFileWriter writer)
    {
        mCapacityStart = capacity;
        setCapacity(capacity);
        mPositionBuffer = positionBuffer;

        mRecords = Lists.newArrayListWithCapacity(mCapacity);
        mSortedRecords = Lists.newArrayListWithCapacity(mCapacity);
        mWriter = writer;
        mCurrentChromosome = "";
        mLastWrittenPosition = 0;
        mMinCachedPosition = 0;
        mReadsWritten = 0;
        mWriteCount = 0;
        mMaxCached = 0;
    }

    public void initialiseStartPosition(final String chromosome, int startPosition)
    {
        mCurrentChromosome = chromosome;
        mLastWrittenPosition = startPosition;
        mUpperBoundPosition = startPosition + mPositionBuffer;
    }

    public void setUpperBoundPosition(int position)
    {
        mUpperBoundPosition = max(position, mLastWrittenPosition + mPositionBuffer);
    }

    private static final double CAPACITY_CHECK_PERCENT = 0.01;
    private static final double CAPACITY_GROW_PERCENT = 0.9;
    private static final double CAPACITY_SHRINK_PERCENT = 0.5;

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

    public void addRecord(final SAMRecord read)
    {
        mRecords.add(read);

        if(mMinCachedPosition == 0 || mRecords.size() < mCapacityThresholdCheck)
        {
            if(mMinCachedPosition == 0)
                mMinCachedPosition = read.getAlignmentStart();
            else
                mMinCachedPosition = min(mMinCachedPosition, read.getAlignmentStart());

            return;
        }

        int maxPosition = mUpperBoundPosition - mPositionBuffer;

        if(maxPosition <= mMinCachedPosition)
            return;

        mMaxCached = max(mMaxCached, mRecords.size());

        // gather all reads prior to this position
        int index = 0;
        mMinCachedPosition = maxPosition;
        while(index < mRecords.size())
        {
            SAMRecord cachedRead = mRecords.get(index);

            if(cachedRead.getAlignmentStart() <= maxPosition)
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

        checkCapacity(maxPosition);
    }

    private void checkCapacity(int maxPosition)
    {
        // check whether capacity needs to grow or shrink
        if(mRecords.size() > mCapacityGrowThreshold)
        {
            MD_LOGGER.debug("sorted read cache pos({}:{}) growing capacity({} -> {})",
                    mCurrentChromosome, maxPosition, mCapacity, mCapacity + mCapacityStart);

            setCapacity(mCapacity + mCapacityStart);
        }
        else if(mCapacity > mCapacityStart && mRecords.size() < mCapacityShrinkThreshold)
        {
            MD_LOGGER.debug("sorted read cache pos({}:{}) reducing capacity({} -> {})",
                    mCurrentChromosome, maxPosition, mCapacity, mCapacityStart);

            setCapacity(mCapacityStart);
        }
    }

    private void setCapacity(int capacity)
    {
        mCapacity = capacity;
        mCapacityGrowThreshold = (int)(capacity * CAPACITY_GROW_PERCENT);
        mCapacityShrinkThreshold = (int)(max(capacity * CAPACITY_SHRINK_PERCENT, 1));
        mCapacityThresholdCheck = (int)max(capacity * CAPACITY_CHECK_PERCENT, 1);
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

    public int averageWriteCount() { return mWriteCount > 0 ? (int)(mReadsWritten / (double)mWriteCount) : 0; }

    public String toString()
    {
        return format("chr(%s) bounds(lastWrite=%d upper=%d) cached(%d/%d max=%d) avgWriteCount(%d from %d) thresholds(check=%d grow=%d shrink=%d)",
                mCurrentChromosome, mLastWrittenPosition, mUpperBoundPosition, mRecords.size(), mCapacity, mMaxCached,
                averageWriteCount(), mWriteCount, mCapacityThresholdCheck, mCapacityGrowThreshold, mCapacityShrinkThreshold);
    }
}
