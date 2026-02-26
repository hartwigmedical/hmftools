package com.hartwig.hmftools.qsee.status;

public enum QcStatusType
{
    NONE,

    PASS,
    WARN,
    FAIL;

    public String displayString()
    {
        return this == QcStatusType.NONE ? "" : this.toString();
    }
}
