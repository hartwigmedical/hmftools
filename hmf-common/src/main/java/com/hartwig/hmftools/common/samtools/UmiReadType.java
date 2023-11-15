package com.hartwig.hmftools.common.samtools;

public enum UmiReadType
{
    NONE,
    SINGLE,
    DUAL;

    @Deprecated
    public static final String DUAL_STRAND_OLD = "DUAL_STRAND";
}
