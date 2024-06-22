package com.hartwig.hmftools.redux.merge;

import static java.lang.String.format;

import java.io.File;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamSequenceReader
{
    private final SamReader mSamReader;
    private final SAMRecordIterator mSamIterator;
    private final String mFilename;

    private SAMRecord mCurrentRecord;
    private int mCurentChromosomeRank;

    public BamSequenceReader(final String refGenomeFile, final String bamFile, final SequenceInfo sequenceInfo)
    {
        File file = new File(bamFile);
        mFilename = file.getName();
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile)).open(file);
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
