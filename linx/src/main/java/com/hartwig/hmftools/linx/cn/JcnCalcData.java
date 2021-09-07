package com.hartwig.hmftools.linx.cn;

public class JcnCalcData
{
    public final double JcnEstimate;
    public final double JcnUncertainty;
    public final boolean Valid;

    public JcnCalcData(double estimate, double uncertainty)
    {
        JcnEstimate = estimate;
        JcnUncertainty = uncertainty;
        Valid = true;
    }

    public JcnCalcData(double estimate, double uncertainty, boolean valid)
    {
        JcnEstimate = estimate;
        JcnUncertainty = uncertainty;
        Valid = valid;
    }
}
