package com.hartwig.hmftools.isofox.common;

public enum GeneMatchType
{
    TOTAL,
    TRANS_SUPPORTING,
    ALT,
    UNSPLICED,
    READ_THROUGH,
    CHIMERIC,
    DUPLICATE,
    MAX;

    public static int typeAsInt(GeneMatchType type)
    {
        switch(type)
        {
            case TOTAL: return 0;
            case TRANS_SUPPORTING: return 1;
            case ALT: return 2;
            case UNSPLICED: return 3;
            case READ_THROUGH: return 4;
            case CHIMERIC: return 5;
            case DUPLICATE: return 6;
            case MAX: return 7;

            default: return 0;
        }
    }
}
