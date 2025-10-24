package com.hartwig.hmftools.common.metrics;

import static java.lang.String.format;

public class ValueFrequency
{
    public final int Value;
    public final long Count;

    public ValueFrequency(final int value, final long count)
    {
        Value = value;
        Count = count;
    }

    public String toString()
    {
        return format("%d=%d", Value, Count);
    }
}
