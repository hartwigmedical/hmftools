package com.hartwig.hmftools.sage.sync;

import static java.lang.Math.min;
import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.S;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class TruncatedBases
{
    public int StartOffset; // bases to be ignored at the start of the combined read
    public int EndOffset;
    public int ForwardStrandStartPosition; // combined fragment start position, taken from the 5' read start

    public static final TruncatedBases NO_TRUNCATION = new TruncatedBases();

    public TruncatedBases()
    {
        StartOffset = 0;
        EndOffset = 0;
        ForwardStrandStartPosition = 0;
    }

    public static TruncatedBases checkShortFragmentTruncation(
            final SAMRecord first, final SAMRecord second, final int firstEffectivePosStart, final int secondEffectivePosStart,
            final int firstEffectivePosEnd, final int secondEffectivePosEnd)
    {
        if(first.getReadNegativeStrandFlag() == second.getReadNegativeStrandFlag())
            return NO_TRUNCATION;

        if(!first.getReadNegativeStrandFlag())
        {
            if(secondEffectivePosStart >= firstEffectivePosStart && firstEffectivePosEnd <= secondEffectivePosEnd)
                return NO_TRUNCATION;
        }
        else
        {
            if(firstEffectivePosStart >= secondEffectivePosStart && secondEffectivePosEnd <= firstEffectivePosEnd)
                return NO_TRUNCATION;
        }

        // truncate bases on the 3' read end which extend past the 5' read start
        TruncatedBases trucatedBases = new TruncatedBases();

        // truncate any bases extending past the 5' read end on the 3' read end
        if(!first.getReadNegativeStrandFlag())
        {
            trucatedBases.ForwardStrandStartPosition = first.getAlignmentStart();

            if(secondEffectivePosStart < firstEffectivePosStart)
            {
                trucatedBases.StartOffset = calculateStartBaseOffset(second, secondEffectivePosStart, firstEffectivePosStart);
                // trucatedBases.StartOffset = firstEffectivePosStart - secondEffectivePosStart;
            }

            if(firstEffectivePosEnd > secondEffectivePosEnd)
            {
                trucatedBases.EndOffset = calculateEndBaseOffset(first, firstEffectivePosEnd, secondEffectivePosEnd);
                // trucatedBases.EndOffset = firstEffectivePosEnd - secondEffectivePosEnd;
            }
        }
        else
        {
            trucatedBases.ForwardStrandStartPosition = second.getAlignmentStart();

            if(firstEffectivePosStart < secondEffectivePosStart)
            {
                trucatedBases.StartOffset = calculateStartBaseOffset(first, firstEffectivePosStart, secondEffectivePosStart);
                // trucatedBases.StartOffset = secondEffectivePosStart - firstEffectivePosStart;
            }

            if(secondEffectivePosEnd > firstEffectivePosEnd)
            {
                trucatedBases.EndOffset = calculateEndBaseOffset(second, secondEffectivePosEnd, firstEffectivePosEnd);
                // trucatedBases.EndOffset = secondEffectivePosEnd - firstEffectivePosEnd;
            }
        }

        return trucatedBases;
    }

    private static int calculateStartBaseOffset(final SAMRecord record, final int effectiveStart, final int actualStart)
    {
        int positionDiff = actualStart - effectiveStart;
        int baseOffset = 0;

        for(CigarElement element : record.getCigar())
        {
            int minBases = min(element.getLength(), positionDiff);

            if(element.getOperator().consumesReferenceBases() || element.getOperator() == S)
                positionDiff -= minBases;

            if(element.getOperator().consumesReadBases())
                baseOffset += minBases;

            if(positionDiff <= 0)
                break;
        }

        return baseOffset;
    }

    private static int calculateEndBaseOffset(final SAMRecord record, final int effectiveEnd, final int actualEnd)
    {
        int positionDiff = effectiveEnd - actualEnd;
        int baseOffset = 0;

        for(int i = record.getCigar().getCigarElements().size() - 1; i >= 0; --i)
        {
            CigarElement element = record.getCigar().getCigarElements().get(i);

            int minBases = min(element.getLength(), positionDiff);

            // soft-clips were factored into the effective position
            if(element.getOperator().consumesReferenceBases() || element.getOperator() == S)
                positionDiff -= minBases;

            if(element.getOperator().consumesReadBases())
                baseOffset += minBases;

            if(positionDiff <= 0)
                break;
        }

        return baseOffset;
    }

    public String toString() { return format("offsets(%d - %d) forwardStart(%d)", StartOffset, EndOffset, ForwardStrandStartPosition); }
}
