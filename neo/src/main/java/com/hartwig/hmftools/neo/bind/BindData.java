package com.hartwig.hmftools.neo.bind;

public class BindData
{
    public final String Allele;
    public final String Peptide;
    public final double Affinity;
    public final String OtherInfo;

    public static final String DELIM = ",";

    public BindData(final String allele, String peptide, double affinity, final String otherInfo)
    {
        Allele = allele;
        Peptide = peptide;
        Affinity = affinity;
        OtherInfo = otherInfo;
    }

    public String toString()
    {
        return String.format("allele(%s) pep(%s) affinity(%.1f) otherInfo(%s)",
                Allele, Peptide, Affinity, OtherInfo);
    }

    public static BindData fromCsv(
            final String data, int alleleIndex, int peptideIndex, int affinityIndex, int otherInfoIndex)
    {
        final String[] items = data.split(DELIM, -1);

        String allele = items[alleleIndex].replaceAll("HLA-", "");
        allele = allele.replaceAll(":", "").replaceAll("\\*", "");

        return new BindData(allele, items[peptideIndex], Double.parseDouble(items[affinityIndex]), items[otherInfoIndex]);
    }
}
