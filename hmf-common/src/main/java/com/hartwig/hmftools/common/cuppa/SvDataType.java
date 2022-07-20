package com.hartwig.hmftools.common.cuppa;

public enum SvDataType
{
    LINE,
    SIMPLE_DEL_20KB_1MB,
    SIMPLE_DUP_32B_200B,
    SIMPLE_DUP_100KB_5MB,
    MAX_COMPLEX_SIZE,
    TELOMERIC_SGL;

    public static int typeIndex(final SvDataType type)
    {
        switch(type)
        {
            case LINE: return 0;
            case SIMPLE_DEL_20KB_1MB: return 1;
            case SIMPLE_DUP_32B_200B: return 2;
            case SIMPLE_DUP_100KB_5MB: return 3;
            case MAX_COMPLEX_SIZE: return 4;
            case TELOMERIC_SGL: return 5;
            default: return -1;
        }
    }

    public static int count() { return SvDataType.values().length; }
}
