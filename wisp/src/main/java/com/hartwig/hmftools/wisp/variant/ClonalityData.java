package com.hartwig.hmftools.wisp.variant;

public class ClonalityData
{
    public final ClonalityMethod Method;
    public final double Vaf;
    public final double VafLow;
    public final double VafHigh;
    public final int VariantCount;

    public final String ExtraInfo;

    public ClonalityData(
            final ClonalityMethod method, final double vaf, final double vafLow, final double vafHigh, int varCount,
            final String extraInfo)
    {
        Method = method;
        Vaf = vaf;
        VafLow = vafLow;
        VafHigh = vafHigh;
        VariantCount = varCount;
        ExtraInfo = extraInfo;
    }

    public static ClonalityData NO_RESULT = new ClonalityData(
            ClonalityMethod.NONE, 0, 0, 0, 0, "");
}
