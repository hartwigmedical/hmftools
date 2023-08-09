package com.hartwig.hmftools.amber;

public enum BamReadMode
{
    DEFAULT,
    CRAM,
    EXP;

    public static final int BAM_REGION_GROUP_MAX = 150000;
    public static final int CRAM_REGION_GROUP_MAX = 10000;

    public static final int CRAM_MIN_GAP_START = 10000;
    public static final int CRAM_MIN_GAP_INCREMENT = 1000;

    public static final int BAM_MIN_GAP_START = 2000;
    public static final int BAM_MIN_GAP_INCREMENT = 200;

    public static int groupMax(final BamReadMode mode)
    {
        switch(mode)
        {
            case CRAM: return CRAM_REGION_GROUP_MAX;
            default: return BAM_REGION_GROUP_MAX;
        }
    }

    public static int minGapStart(final BamReadMode mode)
    {
        switch(mode)
        {
            case CRAM: return CRAM_MIN_GAP_START;
            default: return BAM_MIN_GAP_START;
        }
    }

    public static int minGapIncrement(final BamReadMode mode)
    {
        switch(mode)
        {
            case CRAM: return CRAM_MIN_GAP_INCREMENT;
            default: return BAM_MIN_GAP_INCREMENT;
        }
    }
}
