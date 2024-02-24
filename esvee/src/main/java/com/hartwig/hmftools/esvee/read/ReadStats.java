package com.hartwig.hmftools.esvee.read;

import static java.lang.String.format;

public class ReadStats
{
    public int TotalReads;
    public int FilteredBaseQual;
    public int FilteredMapQual;
    public int PolyGTrimmed;
    public int IndelSoftClipConverted;

    public ReadStats()
    {
        TotalReads = 0;
        FilteredBaseQual = 0;
        FilteredMapQual = 0;
        PolyGTrimmed = 0;
        IndelSoftClipConverted = 0;

    }

    public void merge(final ReadStats other)
    {
        TotalReads += other.TotalReads;
        FilteredBaseQual += other.FilteredBaseQual;
        FilteredMapQual += other.FilteredMapQual;
        PolyGTrimmed += other.PolyGTrimmed;
        IndelSoftClipConverted += other.IndelSoftClipConverted;
    }

    public String toString()
    {
        return format("reads(%d) filter(baseQual=%d mapQual=%d) polyGTrim(%d) indelSoftClip(%d)",
                TotalReads, FilteredBaseQual, FilteredMapQual, PolyGTrimmed, IndelSoftClipConverted);
    }
}
