package com.hartwig.hmftools.redux.write;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.FragmentUtils.readToString;

import java.util.SortedSet;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

public class SortedBamWriter
{
    private final SortedBamConfig mConfig;

    private final int mPositionBuffer;

    private final SortedSet<ReadOrAlignmentStart> mRecords;
    private final SAMFileWriter mWriter;

    // state to control which reads can be processed
    private String mCurrentChromosome; // reads on another chromosome cannot be
    private int mLastWrittenPosition; // last write from the sorted cache - no lower positions can be processed
    private int mUpperBoundPosition; // most recent read, forming the upper bound for accepting reads
    private int mUpperWritablePosition; // upper position on reads ready for writing
    private int mMinCachedPosition;

    // stats
    private int mWriteCount;
    private long mReadCount;
    private long mReadsWritten;
    private int mMaxCached;
    private int mMaxWrite;

    private boolean mPerfDebug;

    private static SAMRecordCoordinateComparator READ_COMPARATOR = new SAMRecordCoordinateComparator();
    // wrap reads to allow for easy splitting SortedSet based on alignment start
    private class ReadOrAlignmentStart implements Comparable<ReadOrAlignmentStart>
    {
        private final long mReadCounter; // so we can store duplicate reads
        public final SAMRecord Read;
        private final int mAlignmentStart;

        public ReadOrAlignmentStart(int alignmentStart)
        {
            mReadCounter = 0;
            Read = null;
            mAlignmentStart = alignmentStart;
        }

        public ReadOrAlignmentStart(long readCount, final SAMRecord read)
        {
            mReadCounter = readCount;
            Read = read;
            mAlignmentStart = 0;
        }

        @Override
        public int compareTo(final ReadOrAlignmentStart other)
        {
            // mAlignmentStart is sorted against Read.getAlignmentStart(). If these are equal then the mAlignmentStart instances are sorted
            // after the reads with the same alignment start.
            if(Read == null && other.Read == null)
            {
                if(mAlignmentStart < other.mAlignmentStart)
                    return -1;

                return mAlignmentStart > other.mAlignmentStart ? 1 : 0;
            }

            if(Read == null)
                return mAlignmentStart < other.Read.getAlignmentStart() ? -1 : 1;

            if(other.Read == null)
                return Read.getAlignmentStart() <= other.mAlignmentStart ? -1 : 1;

            int comparison = READ_COMPARATOR.compare(Read, other.Read);
            if(comparison != 0)
                return comparison;

            if(mReadCounter < other.mReadCounter)
                return -1;

            return mReadCounter > other.mReadCounter ? 1 : 0;
        }
    }

    public SortedBamWriter(final SortedBamConfig config, @Nullable final SAMFileWriter writer)
    {
        mConfig = config;
        mPositionBuffer = mConfig.PositionBuffer;

        mRecords = Sets.newTreeSet();
        mWriter = writer;
        mCurrentChromosome = "";
        mLastWrittenPosition = 0;
        mUpperWritablePosition = 0;
        mMinCachedPosition = 0;
        mReadsWritten = 0;
        mWriteCount = 0;
        mReadCount = 0;
        mMaxCached = 0;
        mMaxWrite = 0;
        mPerfDebug = false;
    }

    public void initialiseStartPosition(final String chromosome, int startPosition)
    {
        if(!mCurrentChromosome.equals(chromosome))
        {
            mCurrentChromosome = chromosome;
            mLastWrittenPosition = 0;
        }

        mMinCachedPosition = startPosition;
        mUpperBoundPosition = startPosition + mPositionBuffer;
        mUpperWritablePosition = startPosition;
    }

    public void setUpperBoundPosition(int position) { mUpperBoundPosition = position; }

    public void setUpperWritablePosition(int position)
    {
        mUpperWritablePosition = position - mConfig.ReadPosCacheBuffer;
    }

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
        mRecords.add(new ReadOrAlignmentStart(mReadCount++, read));
        mMaxCached = max(mMaxCached, mRecords.size());

        if(mPerfDebug && (mReadCount % LOG_COUNT) == 0)
        {
            RD_LOGGER.debug("sorted read cache chr({}:{}) records({}) avgWrite({})",
                    mCurrentChromosome, mLastWrittenPosition, mRecords.size(), avgWriteCount());
        }

        if(mMinCachedPosition == 0)
            mMinCachedPosition = read.getAlignmentStart();
        else
            mMinCachedPosition = min(mMinCachedPosition, read.getAlignmentStart());

        if(mUpperWritablePosition < mMinCachedPosition)
            return;

        // see how many records could be written
        SortedSet<ReadOrAlignmentStart> sortedWritableRecords = mRecords.headSet(new ReadOrAlignmentStart(mUpperWritablePosition));
        if(sortedWritableRecords.size() < mConfig.MinWriteCount)
            return;

        writeRecords(sortedWritableRecords);

        mMinCachedPosition = mRecords.isEmpty() ? mUpperWritablePosition + 1 : mRecords.first().Read.getAlignmentStart();
    }

    private void writeRecords(final SortedSet<ReadOrAlignmentStart> records)
    {
        if(records.isEmpty())
            return;

        mReadsWritten += records.size();

        if(mWriter != null)
        {
            // check the most recent write position against this new batch
            if(records.first().Read.getAlignmentStart() < mLastWrittenPosition)
            {
                SAMRecord nextRecord = records.first().Read;

                RD_LOGGER.error("sorted BAM cache({}) writing earlier(readStart={} vs last={}) from {} records, read: {}",
                        toString(), nextRecord.getAlignmentStart(), mLastWrittenPosition, records.size(), readToString(nextRecord));

                System.exit(1);
            }

            for(ReadOrAlignmentStart readOrAlignmentStart : records)
            {
                mWriter.addAlignment(readOrAlignmentStart.Read);
            }
        }

        mMaxWrite = max(mMaxWrite, records.size());

        mLastWrittenPosition = records.last().Read.getAlignmentStart();

        records.clear();
        ++mWriteCount;
    }

    public void flush()
    {
        // write all cached records
        writeRecords(mRecords);
    }

    public long written() { return mReadsWritten; }
    public long writeCount() { return mWriteCount; }
    public int cached() { return mRecords.size(); }
    public int maxCache() { return mMaxCached; }
    public int maxWrite() { return mMaxWrite; }

    public int avgWriteCount() { return mWriteCount > 0 ? (int)(mReadsWritten / (double)mWriteCount) : 0; }

    public String toString()
    {
        return format("chr(%s) lastWritePos(%d) upper(sort=%d bound=%d) cached(%d max=%d) avgWriteCount(%d from %d)",
                mCurrentChromosome, mLastWrittenPosition, mUpperWritablePosition, mUpperBoundPosition, mRecords.size(), mMaxCached,
                avgWriteCount(), mWriteCount);
    }
}
