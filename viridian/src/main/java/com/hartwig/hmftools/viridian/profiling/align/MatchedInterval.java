package com.hartwig.hmftools.viridian.profiling.align;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// A contiguous run of reference bases matched in an alignment.
public record MatchedInterval(
        // 1-based inclusive.
        int referenceStart,
        int length
)
{
    public MatchedInterval
    {
        if(referenceStart < 1)
        {
            throw new IllegalArgumentException("Invalid reference start: " + referenceStart);
        }
        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid interval length: " + length);
        }
    }

    public static List<MatchedInterval> from(int alignmentStart, Cigar cigar, int[] mismatchPositions)
    {
        List<MatchedInterval> intervals = new ArrayList<>(cigar.numCigarElements() + mismatchPositions.length);
        int blockStart = alignmentStart;
        int nextMismatch = 0;
        for(CigarElement element : cigar.getCigarElements())
        {
            CigarOperator operator = element.getOperator();
            if(operator.consumesReadBases() && operator.consumesReferenceBases())
            {
                int start = blockStart;
                int end = blockStart + element.getLength();
                for(; nextMismatch < mismatchPositions.length && mismatchPositions[nextMismatch] < end; ++nextMismatch)
                {
                    int mismatch = mismatchPositions[nextMismatch];
                    if(mismatch < start)
                    {
                        throw inconsistentMismatches(cigar, mismatchPositions);
                    }
                    if(mismatch > start)
                    {
                        intervals.add(new MatchedInterval(start, mismatch - start));
                    }
                    start = mismatch + 1;
                }
                if(end > start)
                {
                    intervals.add(new MatchedInterval(start, end - start));
                }
            }
            if(operator.consumesReferenceBases())
            {
                blockStart += element.getLength();
            }
        }
        if(nextMismatch < mismatchPositions.length)
        {
            throw inconsistentMismatches(cigar, mismatchPositions);
        }
        return intervals;
    }

    private static IllegalArgumentException inconsistentMismatches(Cigar cigar, int[] mismatchPositions)
    {
        return new IllegalArgumentException(String.format(
                "Mismatch positions %s outside the aligned bases of CIGAR %s",
                Arrays.toString(mismatchPositions), cigar));
    }
}
