package com.hartwig.hmftools.qsee.feature;

public enum NumberFormat
{
    NONE,

    NUMBER,
    PERCENT,
    LOG10;

    public String displayString()
    {
        return (this == NumberFormat.NONE) ? "" : this.toString();
    }
}
