package com.hartwig.hmftools.ctdna.purity.variant;

public class ClonalityResult
{
    public final ClonalityMethod Method;
    public final double Vaf;
    public final double VafLow;
    public final double VafHigh;
    public final int VariantCount;

    public ClonalityResult(final ClonalityMethod method, final double vaf, final double vafLow, final double vafHigh, int varCount)
    {
        Method = method;
        Vaf = vaf;
        VafLow = vafLow;
        VafHigh = vafHigh;
        VariantCount = varCount;
    }

    public static ClonalityResult INVALID_RESULT = new ClonalityResult(ClonalityMethod.NONE, 0, 0, 0, 0);
}
