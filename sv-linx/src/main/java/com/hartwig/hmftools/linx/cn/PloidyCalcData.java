package com.hartwig.hmftools.linx.cn;

public class PloidyCalcData
{
    public final double PloidyEstimate;
    public final double PloidyUncertainty;
    public final boolean Valid;

    public PloidyCalcData(double estimate, double uncertainty, boolean valid)
    {
        PloidyEstimate = estimate;
        PloidyUncertainty = uncertainty;
        Valid = valid;
    }

}
