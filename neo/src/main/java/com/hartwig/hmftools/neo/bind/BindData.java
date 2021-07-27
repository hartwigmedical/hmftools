package com.hartwig.hmftools.neo.bind;

public class BindData
{
    public final String Allele;
    public final String Peptide;
    public final double Affinity;
    public final double PredictedAffinity;
    public final String Source;

    public double Score;
    public double RankPerc;

    public static final String DELIM = ",";
    public static final String RANDOM_SOURCE = "Random";

    public BindData(final String allele, String peptide, double affinity, double predictedAffinity, final String source)
    {
        Allele = allele;
        Peptide = peptide;
        Affinity = affinity;
        PredictedAffinity = predictedAffinity;
        Source = source;

        Score = 0;
        RankPerc = 0;
    }

    public int peptideLength() { return Peptide.length(); }
    public boolean isRandom() { return Source.equals(RANDOM_SOURCE); }
    public boolean isTraining() { return !isRandom(); }

    public String toString()
    {
        return String.format("allele(%s) pep(%s) affinity(%.1f pred=%.1f) source(%s)",
                Allele, Peptide, Affinity, PredictedAffinity, Source);
    }

    public static BindData fromCsv(
            final String data, int alleleIndex, int peptideIndex, int affinityIndex, int predictedIndex, int otherInfoIndex)
    {
        final String[] items = data.split(DELIM, -1);

        String allele = items[alleleIndex].replaceAll("HLA-", "");
        allele = allele.replaceAll(":", "").replaceAll("\\*", "");

        return new BindData(
                allele, items[peptideIndex], Double.parseDouble(items[affinityIndex]),
                Double.parseDouble(items[predictedIndex]), items[otherInfoIndex]);
    }

    public static BindData fromCsv(
            final String data, int alleleIndex, int peptideIndex, int predAffinityIndex, double maxAffinity)
    {
        final String[] items = data.split(DELIM, -1);
        double predictedAffinity = Double.parseDouble(items[predAffinityIndex]);
        return new BindData(items[alleleIndex], items[peptideIndex], maxAffinity, predictedAffinity, RANDOM_SOURCE);
    }

}
