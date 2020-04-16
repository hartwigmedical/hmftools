package com.hartwig.hmftools.isofox.common;

public enum RegionMatchType
{
    NONE,
    EXON_BOUNDARY,  // read matches one exon boundary
    WITHIN_EXON,    // read fully contained within the exon
    EXON_MATCH,     // read fully contained within the exon
    EXON_INTRON;      // reads spanning to unmapped regions where adjacent regions exist


    public static int matchRank(RegionMatchType type)
    {
        switch(type)
        {
            case EXON_INTRON:
                return 1;
            case WITHIN_EXON:
                return 2;
            case EXON_BOUNDARY:
                return 3;
            case EXON_MATCH:
                return 3;
            case NONE:
            default:
                return 0;
        }
    }
}
