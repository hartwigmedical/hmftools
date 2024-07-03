package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.String.format;

public class SomaticPeak
{
    public final double AlleleFrequency;
    public final int Count;

    public SomaticPeak(final double alleleFrequency, final int count)
    {
        AlleleFrequency = alleleFrequency;
        Count = count;
    }

    public String toString()
    {
        return format("vaf(%.3f) count(%d)", AlleleFrequency, Count);
    }
}