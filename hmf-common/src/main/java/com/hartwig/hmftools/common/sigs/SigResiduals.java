package com.hartwig.hmftools.common.sigs;

public class SigResiduals
{
    public double Total;
    public double Percent;
    public double Excess;

    public static final String SIG_MISALLOCATED = "MISALLOC";
    public static final String SIG_UNALLOCATED = "UNALLOC";
    public static final String SIG_EXCESS = "EXCESS";

    public SigResiduals(final double total, final double percent, final double excess)
    {
        Total = total;
        Percent = percent;
        Excess = excess;
    }

    public SigResiduals()
    {
        Total = 0;
        Percent = 0;
        Excess = 0;
    }

    public double unallocated() { return Total - Excess; }

}
