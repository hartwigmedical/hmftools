package com.hartwig.hmftools.bamtools.metrics;

public class Statistics
{
    public final double Mean;
    public final double Median;
    public final int MedianAbsoluteDeviation;
    public final double StandardDeviation;

    public Statistics(final double mean, final double median, final int mad, final double standardDeviation)
    {
        Mean = mean;
        Median = median;
        MedianAbsoluteDeviation = mad;
        StandardDeviation = standardDeviation;
    }
}
