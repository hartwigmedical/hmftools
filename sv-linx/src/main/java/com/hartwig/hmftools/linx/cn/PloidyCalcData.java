package com.hartwig.hmftools.linx.cn;

import com.hartwig.hmftools.linx.types.SvVarData;

public class PloidyCalcData
{
    public final double PloidyEstimate;
    public final double PloidyUncertainty;
    public final boolean Valid;

    public PloidyCalcData(double estimate, double uncertainty)
    {
        PloidyEstimate = estimate;
        PloidyUncertainty = uncertainty;
        Valid = true;
    }

    public PloidyCalcData(double estimate, double uncertainty, boolean valid)
    {
        PloidyEstimate = estimate;
        PloidyUncertainty = uncertainty;
        Valid = valid;
    }

    public PloidyCalcData(final SvVarData var)
    {
        PloidyEstimate = var.ploidy();
        PloidyUncertainty = var.ploidyUncertainty();
        Valid = true;
    }
}
