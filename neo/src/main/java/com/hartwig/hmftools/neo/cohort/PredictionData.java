package com.hartwig.hmftools.neo.cohort;

public class PredictionData
{
    public final String Allele;
    public final int NeId;
    public final String Peptide;
    public final double Affinity;
    public final double Presentation;

    public static final String DELIM = ",";

    public PredictionData(final String allele, int neId, String peptide, double affinity, double presentation)
    {
        Allele = allele;
        NeId = neId;
        Peptide = peptide;
        Affinity = affinity;
        Presentation = presentation;
    }

    public static PredictionData fromCsv(
            final String data, int alleleIndex, int neIdIndex, int peptideIndex, int affinityIndex, int presentationIndex)
    {
        final String[] items = data.split(DELIM, -1);

        String allele = items[alleleIndex].replaceAll("HLA-", "");

        return new PredictionData(
                allele, Integer.parseInt(items[neIdIndex]), items[peptideIndex],
                Double.parseDouble(items[affinityIndex]), Double.parseDouble(items[presentationIndex]));
    }
}
