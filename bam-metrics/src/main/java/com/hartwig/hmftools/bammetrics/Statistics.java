package com.hartwig.hmftools.bammetrics;

public class Statistics
{
    public final double Median;
    public final double StandardDeviation;
    public final double MedianAbsoluteDeviation;

    public Statistics(final double median, final double standardDeviation, final double medianAbsoluteDeviation)
    {
        Median = median;
        StandardDeviation = standardDeviation;
        MedianAbsoluteDeviation = medianAbsoluteDeviation;
    }
}
