package com.hartwig.hmftools.neo.bind;

public class RandomPeptideData
{
    public final String Peptide;
    public final String UpFlank;
    public final String DownFlank;

    public RandomPeptideData(final String peptide, final String upFlank, final String downFlank)
    {
        Peptide = peptide;
        UpFlank = upFlank;
        DownFlank = downFlank;
    }
}
