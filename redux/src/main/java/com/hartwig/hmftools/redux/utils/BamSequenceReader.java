package com.hartwig.hmftools.redux.utils;

import static java.lang.String.format;

import java.io.File;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamSequenceReader
{
    private final SamReader mSamReader;
    private final SAMRecordIterator mSamIterator;
    private final String mFilename;
    private final String mChromosome;

    private SAMRecord mCurrentRecord;
    private boolean mOnUnmappedRecords;

    public BamSequenceReader(final String refGenomeFile, final String bamFile, final SAMSequenceRecord sequence)
    {
        File file = new File(bamFile);
        mFilename = file.getName();
        mChromosome = sequence.getContig();
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile)).open(file);
        mCurrentRecord = null;
        mOnUnmappedRecords = false;

        QueryInterval[] queryIntervals = new QueryInterval[] {
                new QueryInterval(sequence.getSequenceIndex(), 1, sequence.getEnd()) };

        mSamIterator = mSamReader.queryOverlapping(queryIntervals);

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

    public boolean isLowerOrEqualWith(final BamSequenceReader other)
    {
        return !isHigherThan(other);
    }

    public boolean isHigherThan(final BamSequenceReader other)
    {
        if(finished())
        {
            return false;
        }

        if(mOnUnmappedRecords != other.onUnmappedRecords())
        {
            return mOnUnmappedRecords;
        }

        return currentPosition() > other.currentPosition();
    }

    public SAMRecord moveNext()
    {
        if(mSamIterator.hasNext())
        {
            mCurrentRecord = mSamIterator.next();

            if(mCurrentRecord.getReadUnmappedFlag() && mCurrentRecord.getMateUnmappedFlag())
            {
                mOnUnmappedRecords = true;
            }
        }
        else
        {
            mCurrentRecord = null;
        }

        return mCurrentRecord;
    }

    public boolean finished()
    {
        return mCurrentRecord == null;
    }

    public boolean onUnmappedRecords()
    {
        return mOnUnmappedRecords;
    }

    public String toString()
    {
        String state;

        if(finished())
        {
            state = "finished";
        }
        else if(mOnUnmappedRecords)
        {
            state = "unmapped";
        }
        else
        {
            state = format("chr(%s:%d)", mChromosome, currentPosition());
        }

        return format("%s: %s", mFilename, state);
    }
}
