package com.hartwig.hmftools.linx.types;

public class ClusterMetrics
{
    public int DBCount;
    public int ShortDBCount;
    public int ClusterDBCount;
    public long TotalDBLength; // sum of deletion bridges
    public long ChainedLength;

    public long TotalRange; // chained distance without double-count and less any isolated short TIs
    public long TraversedRange; // bases covered by the cluster using ploidy data
    public long TotalDeleted; // bases with traversed range but without cluster CN support
    public double ValidAllelePloidySegmentPerc;

    public ClusterMetrics()
    {
        DBCount = 0;
        ShortDBCount = 0;
        ClusterDBCount = 0;
        TotalDBLength = 0;
        ChainedLength = 0;
        TotalRange = 0;
        TraversedRange = 0;
        TotalDeleted = 0;
        ValidAllelePloidySegmentPerc = 1;
    }

}
