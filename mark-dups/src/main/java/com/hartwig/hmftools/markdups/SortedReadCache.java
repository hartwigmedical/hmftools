package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

public class SortedReadCache
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

    // stats
    private int mWriteCount;
    private long mReadsWritten;

    public SortedReadCache(final int capacity, final int positionBuffer, @Nullable final SAMFileWriter writer)
    {
        mCapacityStart = capacity;
        setCapacity(capacity);
        mPositionBuffer = positionBuffer;

        mRecords = Lists.newArrayListWithCapacity(mCapacity);
        mSortedRecords = Lists.newArrayListWithCapacity(mCapacity);
        mWriter = writer;
        mCurrentChromosome = "";
        mReadsWritten = 0;
        mWriteCount = 0;
    }

    private static final double CAPACITY_CHECK_PERCENT = 0.01;
    private static final double CAPACITY_GROW_PERCENT = 0.9;
    private static final double CAPACITY_SHRINK_PERCENT = 0.5;

    public void addRecord(final SAMRecord read)
    {
        if(!mCurrentChromosome.equals(read.getReferenceName()))
        {
            mCurrentChromosome = read.getReferenceName();
            flush();
        }

        if(mRecords.size() < mCapacityThresholdCheck)
        {
            mRecords.add(read);
            return;
        }

        int maxPosition = read.getAlignmentStart() - mPositionBuffer;

        // gather all reads prior to this position
        int index = 0;
        while(index < mRecords.size())
        {
            SAMRecord cachedRead = mRecords.get(index);

            if(cachedRead.getAlignmentStart() < maxPosition)
            {
                mSortedRecords.add(cachedRead);
                mRecords.remove(index);
            }
            else
            {
                // cannot break since the latter records may have earlier positions
                ++index;
            }
        }

        // sort and write them
        sortAndWrite(mSortedRecords);

        mRecords.add(read);

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
        mCapacityShrinkThreshold = (int)(capacity * CAPACITY_SHRINK_PERCENT);
        mCapacityThresholdCheck = (int)(capacity * CAPACITY_CHECK_PERCENT);
    }

    private void sortAndWrite(final List<SAMRecord> records)
    {
        if(records.isEmpty())
            return;

        mReadsWritten += records.size();
        Collections.sort(records, new SAMRecordCoordinateComparator());

        if(mWriter != null)
            records.forEach(x -> mWriter.addAlignment(x));

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

    public int averageWriteCount() { return mWriteCount > 0 ? (int)(mReadsWritten / (double)mWriteCount) : 0; }

    public String toString()
    {
        return format("chromosome({}) cached({}/{}) avgWriteCount({}) thresholds(check={} grow={} shrink={})",
                mCurrentChromosome, mRecords.size(), mCapacity, averageWriteCount(),
                mCapacityThresholdCheck, mCapacityGrowThreshold, mCapacityShrinkThreshold);
    }
}
