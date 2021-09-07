package com.hartwig.hmftools.neo.epitope;

public class PeptideData
{
    public final String Peptide;
    public final String UpFlank;
    public final String DownFlank;

    public PeptideData(final String peptide, final String upFlank, final String downFlank)
    {
        Peptide = peptide;
        UpFlank = upFlank;
        DownFlank = downFlank;
    }
}
