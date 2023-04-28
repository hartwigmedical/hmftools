package com.hartwig.hmftools.linx.types;

public enum ArmClusterType
{
    UNKNOWN,
    ISOLATED_BE,
    TI_ONLY,
    DSB, // includes simple DELs within the arm-cluster window
    FOLDBACK,
    FOLDBACK_DSB,
    COMPLEX_FOLDBACK,
    COMPLEX_LINE,
    COMPLEX_OTHER,
    SIMPLE_DUP,
    SAME_ORIENT;
}
