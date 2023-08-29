package com.hartwig.hmftools.purple.copynumber;

public class LohCalcData
{
    public long LohBaseCount;
    public long TotalBaseCount; // since covers whole genome, not just a chromosome
    public int Segments;

    public LohCalcData() { this(0, 0, 0); }

    public static final LohCalcData LOH_NONE = new LohCalcData();

    public LohCalcData(final long lohBaseCount, final long totalBaseCount, final int segments)
    {
        LohBaseCount = lohBaseCount;
        TotalBaseCount = totalBaseCount;
        Segments = segments;
    }

    public void add(final LohCalcData other)
    {
        LohBaseCount += other.LohBaseCount;
        TotalBaseCount += other.TotalBaseCount;
        Segments += other.Segments;
    }

    public double lohPercent() { return TotalBaseCount > 0 ? LohBaseCount / (double)TotalBaseCount : 0; }
}
