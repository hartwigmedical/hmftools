package com.hartwig.hmftools.common.cuppa;

public enum SvDataType
{
    LINE,
    SIMPLE_DEL_20KB_1MB,
    SIMPLE_DUP_32B_200B,
    SIMPLE_DUP_100KB_5MB,
    MAX_COMPLEX_SIZE,
    TELOMERIC_SGL;

    public static int typeIndex(final SvDataType type) { return type.ordinal(); }

    public static int count()
    {
        return SvDataType.values().length;
    }
}
