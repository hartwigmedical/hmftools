package com.hartwig.hmftools.wisp.purity.variant;

public class ClonalityData
{
    public final ClonalityMethod Method;
    public final double Vaf;
    public final double VafLow;
    public final double VafHigh;
    public final int VariantCount;
    public final double DropoutRate;

    // currently only for the VAF peak mode
    public double PeakBandwidth;
    public double PeakBandwidthLow;
    public double PeakBandwidthHigh;

    public ClonalityData(
            final ClonalityMethod method, final double vaf, final double vafLow, final double vafHigh, int varCount, double dropoutRate,
            double peakBandwidth, double peakBandwidthLow, double peakBandwidthHigh)
    {
        Method = method;
        Vaf = vaf;
        VafLow = vafLow;
        VafHigh = vafHigh;
        VariantCount = varCount;
        DropoutRate = dropoutRate;
        PeakBandwidth = peakBandwidth;
        PeakBandwidthLow = peakBandwidthLow;
        PeakBandwidthHigh = peakBandwidthHigh;
    }

    public static ClonalityData NO_RESULT = new ClonalityData(
            ClonalityMethod.NONE, 0, 0, 0, 0, 0,
            0, 0, 0);
}
