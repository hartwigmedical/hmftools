package com.hartwig.hmftools.pavereverse.base;

public class ChangeContextData
{
    public final int ExonIndex;
    public final int ChangeStart;
    public final int ChangeEnd;
    public final int PaddingInPreviousExon;
    public final int PaddingInNextExon;
    public final int AminoAcidNumberOfFirstAminoAcidStartingInExon;

    public ChangeContextData(final int exonIndex,
            final int changeStart,
            final int changeEnd,
            final int paddingInPreviousExon,
            final int paddingInNextExon,
            final int numberOfCodonsStartingInPreviousExons)
    {
        ExonIndex = exonIndex;
        ChangeStart = changeStart;
        ChangeEnd = changeEnd;
        PaddingInPreviousExon = paddingInPreviousExon;
        PaddingInNextExon = paddingInNextExon;
        AminoAcidNumberOfFirstAminoAcidStartingInExon = numberOfCodonsStartingInPreviousExons + 1;
    }
}
