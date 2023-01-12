package com.hartwig.hmftools.isofox.common;

public enum RegionMatchType
{
    NONE,
    EXON_BOUNDARY,  // read matches one exon boundary
    WITHIN_EXON,    // read fully contained within the exon, ie not touching either boundary
    EXON_INTRON,    // reads spanning to unmapped regions where adjacent regions exist
    INTRON;

    public static boolean validExonMatch(RegionMatchType type)
    {
        return type == EXON_BOUNDARY || type == WITHIN_EXON;
    }

    public static boolean exonBoundary(RegionMatchType type)
    {
        return type == EXON_BOUNDARY;
    }

    public static int matchRank(RegionMatchType type)
    {
        switch(type)
        {
            case EXON_BOUNDARY:
                return 4;
            case WITHIN_EXON:
                return 3;
            case EXON_INTRON:
                return 2;
            case INTRON:
                return 1;
            case NONE:
            default:
                return 0;
        }
    }
}
