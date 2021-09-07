package com.hartwig.hmftools.common.stats;

public class AucData
{
    public final boolean IsPositive;
    public final double Value;

    public final boolean IsPercentile;

    public AucData(final boolean isPositive, final double value, boolean isPercentile)
    {
        IsPositive = isPositive;
        Value = value;
        IsPercentile = isPercentile;
    }

    public String toString() { return String.format("%s=%g result=%s", IsPercentile ? "percentile" : "value", Value, IsPositive); }
}
