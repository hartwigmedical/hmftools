package com.hartwig.hmftools.common.cobalt;

public class MedianRatio
{
    public final String Chromosome;
    public final double MedianRatio;
    public final int Count;

    public MedianRatio(final String chromosome, final double medianRatio, final int count)
    {
        Chromosome = chromosome;
        MedianRatio = medianRatio;
        Count = count;
    }
}
