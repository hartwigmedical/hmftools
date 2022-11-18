package com.hartwig.hmftools.bamtools.metrics;

public class Statistics
{
    public final double Mean;
    public final double Median;
    public final double StandardDeviation;
    public final double MedianAbsoluteDeviation;

    public Statistics(final double mean, final double median, final double standardDeviation, final double medianAbsoluteDeviation)
    {
        Mean = mean;
        Median = median;
        StandardDeviation = standardDeviation;
        MedianAbsoluteDeviation = medianAbsoluteDeviation;
    }
}
