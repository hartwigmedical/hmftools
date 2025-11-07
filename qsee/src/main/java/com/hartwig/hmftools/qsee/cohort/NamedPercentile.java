package com.hartwig.hmftools.qsee.cohort;

public enum NamedPercentile
{
    PCT_MIN("Min", 5),
    PCT_LOWER("Lower", 25),
    PCT_MID("Mid", 50),
    PCT_UPPER("Upper", 75),
    PCT_MAX("Max", 95);

    private final String mSimpleName;
    private final double mPercentile;

    NamedPercentile(String simpleName, double percentile)
    {
        mSimpleName = simpleName;
        mPercentile = percentile;
    }

    public String simpleName() { return mSimpleName; }
    public double percentile() { return mPercentile; }
}
