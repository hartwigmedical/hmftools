package com.hartwig.hmftools.neo.predict;

public class PredictionValues
{
    public final double Affinity;
    public final double AffinityPercentile;
    public final double ProcessingScore;
    public final double PresentationScore;
    public final double PresentationPercentile;

    public PredictionValues(final double affinity, final double affinityPercentile, final double processingScore,
            final double presentationScore, final double presentationPercentile)
    {
        Affinity = affinity;
        AffinityPercentile = affinityPercentile;
        ProcessingScore = processingScore;
        PresentationScore = presentationScore;
        PresentationPercentile = presentationPercentile;
    }

    public static PredictionValues fromCsv(final String[] items, int startIndex)
    {
        return new PredictionValues(
                Double.parseDouble(items[startIndex++]), Double.parseDouble(items[startIndex++]), Double.parseDouble(items[startIndex++]),
                Double.parseDouble(items[startIndex++]), Double.parseDouble(items[startIndex++]));
    }
}
