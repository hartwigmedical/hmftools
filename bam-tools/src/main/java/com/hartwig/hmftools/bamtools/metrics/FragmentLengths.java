package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class FragmentLengths
{
    private int mMaxReadLength;
    private final List<LengthFrequency> mLengthFrequencies;

    public static final int FRAG_LENGTH_DIST_MAX_LENGTH = 1500;

    public FragmentLengths()
    {
        mLengthFrequencies = Lists.newArrayList();
    }

    public List<LengthFrequency> lengthFrequencies() { return mLengthFrequencies; }

    public void processRead(final SAMRecord read)
    {
        // cull invalid reads without waiting for the paired read
        if(!isCandidateRecord(read))
            return;

        mMaxReadLength = max(mMaxReadLength, read.getReadBases().length);

        addFragmentLength(read);
    }

    private boolean isCandidateRecord(final SAMRecord read)
    {
        if(read.getDuplicateReadFlag())
            return false;

        if(read.isSecondaryOrSupplementary())
            return false;

        if(read.getReadPairedFlag())
        {
            if(read.getSecondOfPairFlag())
                return false;

            // ignore translocations and inversions
            if(!read.getMateReferenceName().equals(read.getReferenceName()))
                return false;

            if(mateNegativeStrand(read) == read.getReadNegativeStrandFlag())
                return false;
        }

        // ignore reads with splits
        if(read.getCigar().containsOperator(CigarOperator.N))
            return false;

        return true;
    }

    private void addFragmentLength(final SAMRecord read)
    {
        int fragmentLength = getLengthBucket(abs(read.getInferredInsertSize()));

        if(fragmentLength <= 0)
            return;

        if(fragmentLength > FRAG_LENGTH_DIST_MAX_LENGTH)
            return;

        addLengthFrequency(fragmentLength, 1);
    }

    private void addLengthFrequency(int fragmentLength, int frequency)
    {
        int index = 0;
        while(index < mLengthFrequencies.size())
        {
            final LengthFrequency fragLengthCount = mLengthFrequencies.get(index);

            if(fragLengthCount.Length < fragmentLength)
            {
                ++index;
                continue;
            }

            if(fragLengthCount.Length == fragmentLength)
            {
                fragLengthCount.Frequency += frequency;
                return;
            }

            break;
        }

        LengthFrequency newFragLengthCount = new LengthFrequency(fragmentLength, frequency);
        mLengthFrequencies.add(index, newFragLengthCount);
    }

    private int getLengthBucket(int fragmentLength)
    {
        // round to nearest unit up to 1000, then 10s up to 3000 then 100s
        if(fragmentLength < 1000)
            return fragmentLength;

        if(fragmentLength < 3000)
            return 10 * (int)round(fragmentLength/10.0);

        return 100 * (int)round(fragmentLength/100.0);
    }

    public void merge(final FragmentLengths other)
    {
        other.lengthFrequencies().forEach(x -> addLengthFrequency(x.Length, x.Frequency));
    }
}
