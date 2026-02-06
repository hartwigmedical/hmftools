package com.hartwig.hmftools.common.cider;

public enum IgTcrRegion
{
    V_REGION,
    D_REGION,
    J_REGION,
    CONSTANT;

    public boolean isVJ()
    {
        return this == V_REGION || this == J_REGION;
    }

    public boolean isVDJ()
    {
        return isVJ() || this == D_REGION;
    }
}
