package com.hartwig.hmftools.neo.utils;

public class BindStudyData
{
    public final String Allele;
    public final String Study;
    public final Character AminoAcid;
    public final int Position;
    public final int BindCount;

    public BindStudyData(final String allele, final String study, final char aminoAcid, final int position, final int bindCount)
    {
        Allele = allele;
        Study = study;
        AminoAcid = aminoAcid;
        Position = position;
        BindCount = bindCount;
    }
}
