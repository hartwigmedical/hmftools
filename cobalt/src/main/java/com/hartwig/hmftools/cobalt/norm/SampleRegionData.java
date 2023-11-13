package com.hartwig.hmftools.cobalt.norm;

import static java.lang.String.format;

public class SampleRegionData
{
    public final double ReadDepth;
    public final double GcRatioPanel;
    public final double GcRatioWgs;

    // calculated values
    private double mAdjustedGcRatio;

    public SampleRegionData(final double readDepth, final double gcRatioPanel, final double gcRatioWgs)
    {
        ReadDepth = readDepth;
        GcRatioPanel = gcRatioPanel;
        GcRatioWgs = gcRatioWgs;
        mAdjustedGcRatio = 0;
    }

    public double adjustedGcRatio() { return mAdjustedGcRatio; }
    public void setAdjustedGcRatio(double adjustedGcRatio) { mAdjustedGcRatio = adjustedGcRatio; }

    public String toString()
    {
        return format("reads(%.3f) ratio(panel=%.3f wgs=%f) adjRatio(%.4f)", ReadDepth, GcRatioPanel, GcRatioWgs, mAdjustedGcRatio);
    }
}
