package com.hartwig.hmftools.pavereverse.variants;

import java.util.List;
import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.base.ChangeContextBuilder;
import com.hartwig.hmftools.pavereverse.base.ChangeContextData;

import org.jetbrains.annotations.NotNull;

public class CodonWindow
{
    private final int Start;
    private final int Length;
    private final int End;

    public CodonWindow(final int startIndex1Based, final int numberOfCodons)
    {
        Preconditions.checkArgument(startIndex1Based >= 1);
        Preconditions.checkArgument(numberOfCodons >= 1);
        this.Start = (startIndex1Based - 1) * 3;
        this.Length = numberOfCodons * 3;
        End = Start + Length - 1;
    }

    @NotNull
    public ChangeContextBuilder seekExonLocation(List<Integer> exonLengths)
    {
        int lengthIncludingCurrent = 0;
        for(int i = 0; i < exonLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            int aminoAcidsStartingInPreviousExons = lengthUpToCurrent / 3;
            final Integer exonLength = exonLengths.get(i);
            lengthIncludingCurrent += exonLength;
            if(Start < lengthIncludingCurrent)
            {
                int paddingInPreviousExon = lengthUpToCurrent % 3;
                if(paddingInPreviousExon != 0)
                {
                    aminoAcidsStartingInPreviousExons += 1;
                }
                // The window starts in this exon.
                // We first need to handle the special case in which the window
                // starts at the very end of this exon. In this situation we return
                // data for a change context based in the next exon, but with
                // a 'companion' change context based in this exon.
                int distanceFromWindowStartToCurrentExonEnd = lengthIncludingCurrent - Start;
                if(distanceFromWindowStartToCurrentExonEnd == 2 || distanceFromWindowStartToCurrentExonEnd == 1)
                {
                    // The window starts at the second last or last base of the current exon
                    // and is mainly in the next exon.
                    if(i == exonLengths.size() - 1)
                    {
                        throw new IllegalArgumentException("No exon covering window: " + this);
                    }
                    int startInThisExon = exonLength - distanceFromWindowStartToCurrentExonEnd;
                    ChangeContextData companionData = new ChangeContextData(
                            i,
                            startInThisExon,
                            exonLength - 1,
                            paddingInPreviousExon,
                            overflowIntoNextExon(lengthIncludingCurrent),
                            aminoAcidsStartingInPreviousExons);

                    int lengthIncludingNextExon = lengthIncludingCurrent + exonLengths.get(i + 1);
                    int paddingInNextExon = overflowIntoNextExon(lengthIncludingNextExon);
                    aminoAcidsStartingInPreviousExons = (lengthIncludingCurrent / 3) + 1;
                    ChangeContextData data = new ChangeContextData(
                            i + 1,
                            0,
                            Length - distanceFromWindowStartToCurrentExonEnd - 1,
                            distanceFromWindowStartToCurrentExonEnd, paddingInNextExon, aminoAcidsStartingInPreviousExons);
                    return new ChangeContextBuilder(data, companionData);
                }
                // The window may extend 1 or 2 bases into the next exon.
                int distanceOfWindowEndBeyondCurrentExon = End - lengthIncludingCurrent + 1;
                if(distanceOfWindowEndBeyondCurrentExon > 2)
                {
                    throw new IllegalArgumentException("Window end overlaps by more than one codon with next exon: " + this);
                }
                if(distanceOfWindowEndBeyondCurrentExon == 1 || distanceOfWindowEndBeyondCurrentExon == 2)
                {
                    return new ChangeContextBuilder(i, Start - lengthUpToCurrent,
                            exonLength - 1, paddingInPreviousExon, distanceOfWindowEndBeyondCurrentExon, aminoAcidsStartingInPreviousExons);
                }
                // The window is entirely within this exon.
                int paddingInNextExon = overflowIntoNextExon(lengthIncludingCurrent);
                int startIndexInExon = Start - lengthUpToCurrent;
                int stopIndex = startIndexInExon + Length - 1;
                return new ChangeContextBuilder(i, startIndexInExon, stopIndex, paddingInPreviousExon, paddingInNextExon, aminoAcidsStartingInPreviousExons);
            }
        }
        throw new IllegalArgumentException("No exon covering window: " + this);
    }

    private static int overflowIntoNextExon(int lengthToEndOfCurrentExon)
    {
        return (3 - lengthToEndOfCurrentExon % 3) % 3;
    }

    @Override
    public String toString()
    {
        return "CodonWindow{" +
                "Start=" + Start +
                ", Length=" + Length +
                '}';
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final CodonWindow that = (CodonWindow) o;
        return Start == that.Start && Length == that.Length;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Start, Length);
    }
}
