package com.hartwig.hmftools.pavereverse.base;

public class ChangeContextData
{
    public final int ExonIndex;
    public final int ChangeStart;
    public final int ChangeEnd;
    public final int PaddingInPreviousExon;
    public final int PaddingInNextExon;
    public final int AminoAcidNumberOfFirstAminoAcidStartingInExon;

    public ChangeContextData(int exonIndex,
            int changeStart,
            int changeEnd,
            int paddingInPreviousExon,
            int paddingInNextExon,
            int numberOfCodonsStartingInPreviousExons)
    {
        ExonIndex = exonIndex;
        ChangeStart = changeStart;
        ChangeEnd = changeEnd;
        PaddingInPreviousExon = paddingInPreviousExon;
        PaddingInNextExon = paddingInNextExon;
        AminoAcidNumberOfFirstAminoAcidStartingInExon = numberOfCodonsStartingInPreviousExons + 1;
    }
}
