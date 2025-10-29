package com.hartwig.hmftools.cobalt.norm;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GCBucket.calcGcBucket;

public class SampleRegionData
{
    public final double ReadDepth;
    public final double PanelGcContent;
    public final int PanelGcBucket;
    public final double PanelGcRatio; // only used when writing temp data to file
    public final double WgsGcRatio;

    // calculated values
    private double mAdjustedGcRatio;

    public SampleRegionData(final double readDepth, final double panelGcContent, final double panelGcRatio, final double wgsGcRatio)
    {
        ReadDepth = readDepth;
        PanelGcRatio = panelGcRatio;
        PanelGcContent = panelGcContent;
        PanelGcBucket = calcGcBucket(PanelGcContent);
        WgsGcRatio = wgsGcRatio;
        mAdjustedGcRatio = 0;
    }

    public double adjustedGcRatio() { return mAdjustedGcRatio; }
    public void setAdjustedGcRatio(double adjustedGcRatio) { mAdjustedGcRatio = adjustedGcRatio; }

    public String toString()
    {
        return format("reads(%.3f) ratio(%.3f wgs=%f) adjRatio(%.4f)", ReadDepth, PanelGcContent, WgsGcRatio, mAdjustedGcRatio);
    }
}
