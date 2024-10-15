package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.String.format;

public class ReadStats
{
    public int TotalReads;
    public int LowBaseQualFiltered;
    public int PolyGTrimmed;
    public int LowBaseQualTrimmed;
    public int IndelSoftClipConverted;
    public int DecoySequences;
    public int RefBaseAlignmentFails;
    public int IdenticalSupplementaries;

    public ReadStats()
    {
        TotalReads = 0;
        LowBaseQualFiltered = 0;
        PolyGTrimmed = 0;
        LowBaseQualTrimmed = 0;
        IndelSoftClipConverted = 0;
        DecoySequences = 0;
        RefBaseAlignmentFails = 0;
        IdenticalSupplementaries = 0;
    }

    public void merge(final ReadStats other)
    {
        TotalReads += other.TotalReads;
        LowBaseQualFiltered += other.LowBaseQualFiltered;
        PolyGTrimmed += other.PolyGTrimmed;
        LowBaseQualTrimmed += other.LowBaseQualTrimmed;
        IndelSoftClipConverted += other.IndelSoftClipConverted;
        DecoySequences += other.DecoySequences;
        RefBaseAlignmentFails += other.RefBaseAlignmentFails;
        IdenticalSupplementaries += other.IdenticalSupplementaries;
    }

    public String toString()
    {
        return format("reads(%d) lowBaseQual(%d) trim(polyG=%d lowBase=%d) indelSoftClip(%d) decoySequences(%d) refBaseAlignFails(%d)",
                TotalReads, LowBaseQualFiltered, PolyGTrimmed, LowBaseQualTrimmed, IndelSoftClipConverted, DecoySequences,
                RefBaseAlignmentFails);
    }
}
