package com.hartwig.hmftools.common.variant.pon;

import static java.lang.String.format;

public class PonFilters
{
    public final int SampleCount;
    public final int MaxReadCount;

    public PonFilters(final int sampleCount, final int maxReadCount)
    {
        SampleCount = sampleCount;
        MaxReadCount = maxReadCount;
    }

    public String toString()
    {
        return format("samples(%d) maxReads(%d)", SampleCount, MaxReadCount);
    }
}
