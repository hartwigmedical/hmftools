package com.hartwig.hmftools.pavereverse.gene;

import java.util.List;

public class ExonForLocation
{
    public final int LengthUpToCurrent;
    public final int LengthIncludingCurrent;
    public final int ExonLength;
    public final int ExonIndex;
    public final int LocationInExon;

    public ExonForLocation(int location, List<Integer> exonLengths)
    {
        int lengthUpToCurrent = 0;
        int lengthIncludingCurrent = 0;
        int exonLength = 0;
        int exonIndex = -1;
        for(int i = 0; i < exonLengths.size(); i++)
        {
            lengthUpToCurrent = lengthIncludingCurrent;
            exonLength = exonLengths.get(i);
            lengthIncludingCurrent += exonLength;
            if (location < lengthIncludingCurrent) {
                exonIndex = i;
                break;
            }
        }
        if(exonIndex < 0)
        {
            throw new IllegalArgumentException("No exon containing positiong: " + location);
        }
        LengthUpToCurrent = lengthUpToCurrent;
        ExonLength = exonLength;
        LengthIncludingCurrent = lengthIncludingCurrent;
        ExonIndex = exonIndex;
        LocationInExon = location - LengthUpToCurrent;
    }

    @Override
    public String toString()
    {
        return "ExonForLocation{" +
                "LengthUpToCurrent=" + LengthUpToCurrent +
                ", LengthIncludingCurrent=" + LengthIncludingCurrent +
                ", ExonLength=" + ExonLength +
                ", ExonIndex=" + ExonIndex +
                '}';
    }
}
