package com.hartwig.hmftools.common.stats;

public class AucData
{
    public final boolean IsPositive;
    public final double Value;

    public AucData(final boolean isPositive, final double value)
    {
        IsPositive = isPositive;
        Value = value;
    }

    public String toString() { return String.format("value=%g result=%s", Value, IsPositive); }
}
