package com.hartwig.hmftools.cobalt.utils;

import static java.lang.String.format;

public class SampleRegionData
{
    public final int ReadCount;
    public final double GcRatioPanel;
    public final double GcRatioWgs;

    // calculated values
    private double mAdjustedGcRatio;

    public SampleRegionData(final int readCount, final double gcRatioPanel, final double gcRatioWgs)
    {
        ReadCount = readCount;
        GcRatioPanel = gcRatioPanel;
        GcRatioWgs = gcRatioWgs;
        mAdjustedGcRatio = 0;
    }

    public double adjustedGcRatio() { return mAdjustedGcRatio; }
    public void setAdjustedGcRatio(double adjustedGcRatio) { mAdjustedGcRatio = adjustedGcRatio; }

    public String toString()
    {
        return format("reads(%d) ratio(panel=%.3f wgs=%f) adjRatio(%.4f)", ReadCount, GcRatioPanel, GcRatioWgs, mAdjustedGcRatio);
    }
}
