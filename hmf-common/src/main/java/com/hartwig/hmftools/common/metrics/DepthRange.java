package com.hartwig.hmftools.common.metrics;

public record DepthRange(int lower, int upper)
{
    public int midpoint()
    {
        return (lower + upper + 1) / 2;
    }
}
