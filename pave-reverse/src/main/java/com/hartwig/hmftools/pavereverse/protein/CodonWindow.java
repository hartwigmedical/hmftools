package com.hartwig.hmftools.pavereverse.protein;

import java.util.List;
import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.base.ChangeContextBuilder;
import com.hartwig.hmftools.pavereverse.base.ChangeContextData;
import com.hartwig.hmftools.pavereverse.gene.ExonForLocation;

public class CodonWindow
{
    private final int mStart;
    private final int mLength;
    private final int mEnd;

    public CodonWindow(int startIndex1Based, int numberOfCodons)
    {
        Preconditions.checkArgument(startIndex1Based >= 1);
        Preconditions.checkArgument(numberOfCodons >= 1);
        mStart = (startIndex1Based - 1) * 3;
        mLength = numberOfCodons * 3;
        mEnd = mStart + mLength - 1;
    }

    public ChangeContextBuilder seekExonLocation(List<Integer> exonLengths)
    {
        ExonForLocation exonForLocation = new ExonForLocation(mStart, exonLengths);
        int aminoAcidsStartingInPreviousExons = exonForLocation.LengthUpToCurrent / 3;
        int paddingInPreviousExon = exonForLocation.LengthUpToCurrent % 3;
        if(paddingInPreviousExon != 0)
        {
            aminoAcidsStartingInPreviousExons += 1;
        }
        // We first need to handle the special case in which the window
        // starts at the very end of this exon. In this situation we return
        // data for a change context based in the next exon, but with
        // a 'companion' change context based in this exon.
        int distanceFromWindowStartToCurrentExonEnd = exonForLocation.LengthIncludingCurrent - mStart;
        if(distanceFromWindowStartToCurrentExonEnd == 2 || distanceFromWindowStartToCurrentExonEnd == 1)
        {
            // The window starts at the second last or last base of the current exon
            // and is mainly in the next exon.
            if(exonForLocation.ExonIndex == exonLengths.size() - 1)
            {
                throw new IllegalArgumentException("No exon covering window: " + this);
            }
            int startInThisExon = exonForLocation.ExonLength - distanceFromWindowStartToCurrentExonEnd;
            ChangeContextData companionData = new ChangeContextData(
                    exonForLocation.ExonIndex,
                    startInThisExon,
                    exonForLocation.ExonLength - 1,
                    paddingInPreviousExon,
                    overflowIntoNextExon(exonForLocation.LengthIncludingCurrent),
                    aminoAcidsStartingInPreviousExons);

            int lengthIncludingNextExon = exonForLocation.LengthIncludingCurrent + exonLengths.get(exonForLocation.ExonIndex + 1);
            int paddingInNextExon = overflowIntoNextExon(lengthIncludingNextExon);
            aminoAcidsStartingInPreviousExons = (exonForLocation.LengthIncludingCurrent / 3) + 1;
            ChangeContextData data = new ChangeContextData(
                    exonForLocation.ExonIndex + 1,
                    0,
                    mLength - distanceFromWindowStartToCurrentExonEnd - 1,
                    distanceFromWindowStartToCurrentExonEnd, paddingInNextExon, aminoAcidsStartingInPreviousExons);
            return new ChangeContextBuilder(data, companionData);
        }
        // The window may extend 1 or 2 bases into the next exon.
        int distanceOfWindowEndBeyondCurrentExon = mEnd - exonForLocation.LengthIncludingCurrent + 1;
        if(distanceOfWindowEndBeyondCurrentExon > 2)
        {
            throw new IllegalArgumentException("Window end overlaps by more than one codon with next exon: " + this);
        }
        if(distanceOfWindowEndBeyondCurrentExon == 1 || distanceOfWindowEndBeyondCurrentExon == 2)
        {
            return new ChangeContextBuilder(exonForLocation.ExonIndex, mStart - exonForLocation.LengthUpToCurrent,
                    exonForLocation.ExonLength
                            - 1, paddingInPreviousExon, distanceOfWindowEndBeyondCurrentExon, aminoAcidsStartingInPreviousExons);
        }
        // The window is entirely within this exon.
        int paddingInNextExon = overflowIntoNextExon(exonForLocation.LengthIncludingCurrent);
        int startIndexInExon = exonForLocation.LocationInExon;
        int stopIndex = startIndexInExon + mLength - 1;
        return new ChangeContextBuilder(exonForLocation.ExonIndex, startIndexInExon, stopIndex, paddingInPreviousExon, paddingInNextExon, aminoAcidsStartingInPreviousExons);
    }

    private static int overflowIntoNextExon(int lengthToEndOfCurrentExon)
    {
        return (3 - lengthToEndOfCurrentExon % 3) % 3;
    }

    @Override
    public String toString()
    {
        return "CodonWindow{" +
                "Start=" + mStart +
                ", Length=" + mLength +
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
        return mStart == that.mStart && mLength == that.mLength;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mStart, mLength);
    }
}
