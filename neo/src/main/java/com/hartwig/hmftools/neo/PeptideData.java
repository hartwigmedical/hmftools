package com.hartwig.hmftools.neo;

public class PeptideData
{
    public final String Peptide;
    public final String UpFlank;
    public final String DownFlank;
    public final double TPM;

    public PeptideData(final String peptide, final String upFlank, final String downFlank)
    {
        this(peptide, upFlank, downFlank, -1);
    }

    public PeptideData(final String peptide, final String upFlank, final String downFlank, final double tpm)
    {
        Peptide = peptide;
        UpFlank = upFlank;
        DownFlank = downFlank;
        TPM = tpm;
    }
}
