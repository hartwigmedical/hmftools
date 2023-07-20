package com.hartwig.hmftools.purple.fitting;

public class SomaticPeak {

    public final double AlleleFrequency;
    public final int Count;

    public SomaticPeak(final double alleleFrequency, final int count)
    {
        AlleleFrequency = alleleFrequency;
        Count = count;
    }
}
