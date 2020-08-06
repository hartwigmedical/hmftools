package com.hartwig.hmftools.cup.svs;

public enum SvDataType
{
    LINE,
    FRAGILE_SITE,
    DUP_SHORT,
    DUP_MEDIUM,
    DUP_LONG,
    MAX_EVENT_SIZE,
    TELOMERIC_SGL;

    public static int typeIndex(final SvDataType type)
    {
        switch(type)
        {
            case LINE: return 0;
            case FRAGILE_SITE: return 1;
            case DUP_SHORT: return 2;
            case DUP_MEDIUM: return 3;
            case DUP_LONG: return 4;
            case MAX_EVENT_SIZE: return 5;
            case TELOMERIC_SGL: return 6;
            default: return -1;
        }
    }

    public static int count() { return SvDataType.values().length; }
}
