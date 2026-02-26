package com.hartwig.hmftools.qsee.feature;

public enum NumberFormat
{
    NONE,

    NUMBER,
    PERCENT,
    LOG;

    public String displayString()
    {
        return (this == NumberFormat.NONE) ? "" : this.toString();
    }
}
