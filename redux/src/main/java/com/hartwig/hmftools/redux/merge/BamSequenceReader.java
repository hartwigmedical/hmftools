package com.hartwig.hmftools.redux.merge;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import static htsjdk.samtools.ValidationStringency.SILENT;

import java.io.File;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamSequenceReader
{
    private final SamReader mSamReader;
    private final SAMRecordIterator mSamIterator;
    private final String mFilename;

    private final SequenceInfo mSequenceInfo;
    private boolean mCheckRecordsStart;

    private SAMRecord mCurrentRecord;
    private int mCurentChromosomeRank;

    public BamSequenceReader(final String refGenomeFile, final String bamFile, final SequenceInfo sequenceInfo)
    {
        File file = new File(bamFile);
        mFilename = file.getName();
        mSamReader = SamReaderFactory.makeDefault()
                .validationStringency(SILENT)
                .referenceSequence(new File(refGenomeFile)).open(file);

        mSequenceInfo = sequenceInfo;

        mCheckRecordsStart = sequenceInfo.Intervals.get(0).start > 1;
        mCurrentRecord = null;
        mCurentChromosomeRank = -1;

        mSamIterator = mSamReader.queryOverlapping(sequenceInfo.asArray());
        moveNext();
    }

    public String filename()
    {
        return mFilename;
    }

    public SAMRecord current()
    {
        return mCurrentRecord;
    }

    public int currentPosition()
    {
        return mCurrentRecord != null ? mCurrentRecord.getAlignmentStart() : -1;
    }
    public int currentChromosomeRank() { return mCurentChromosomeRank; }

    public boolean isLowerOrEqualWith(final BamSequenceReader other) { return !isHigherThan(other); }

    public boolean isHigherThan(final BamSequenceReader other)
    {
        if(finished())
            return false;

        if(mCurentChromosomeRank != other.currentChromosomeRank())
            return mCurentChromosomeRank > other.currentChromosomeRank();
        else
            return currentPosition() > other.currentPosition();
    }

    public SAMRecord moveNext()
    {
        if(mCheckRecordsStart)
        {
            // ignore reads starting before the first interval - all subsequent intervals will start at 1 in the next chromosome
            mCheckRecordsStart = false;

            QueryInterval firstInterval = mSequenceInfo.Intervals.get(0);

            int skippedRecords = 0;

            while(mSamIterator.hasNext())
            {
                mCurrentRecord = mSamIterator.next();

                if(mCurrentRecord.getReferenceIndex() == firstInterval.referenceIndex && mCurrentRecord.getAlignmentStart() < firstInterval.start)
                {
                    ++skippedRecords;
                    continue;
                }

                if(skippedRecords > 0)
                {
                    RD_LOGGER.trace("seqRangeId({}) skipped {} records to start of first interval({})",
                            mSequenceInfo.Id, skippedRecords, firstInterval);
                }

                mCurentChromosomeRank = mCurrentRecord.getReferenceIndex();
                return mCurrentRecord;
            }
        }

        if(mSamIterator.hasNext())
        {
            mCurrentRecord = mSamIterator.next();

            if(mCurentChromosomeRank != mCurrentRecord.getReferenceIndex())
            {
                mCurentChromosomeRank = mCurrentRecord.getReferenceIndex();
            }
        }
        else
        {
            mCurrentRecord = null;
            mCurentChromosomeRank = -1;
        }

        return mCurrentRecord;
    }

    public boolean finished()
    {
        return mCurrentRecord == null;
    }

    public String toString()
    {
        String state;

        if(finished())
        {
            state = "finished";
        }
        else
        {
            state = format("chr(%s:%d)", mCurrentRecord.getReferenceName(), currentPosition());
        }

        return format("%s: %s", mFilename, state);
    }
}
