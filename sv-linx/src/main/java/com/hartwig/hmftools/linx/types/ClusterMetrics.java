package com.hartwig.hmftools.linx.types;

public class ClusterMetrics
{
    public int DBCount;
    public int ShortDBCount;
    public int ClusterDBCount;

    // sum of deletion bridges including regions interupted by other SVs (eg unphased)
    public long TotalDBLength;

    // where (typically long) TIs cross unclustered deleted regions (typically simple DELs)
    // an indication that the cluster could have included these in shattering events
    public int TraversedDelLength;
    public int TraversedDelCount;

    public long ChainedLength; // sum of chain lengths, allowing for traversal of the same genomic region (ie from a DUP)

    public long TotalRange; // chained distance without double-count and less any isolated short TIs
    public long TraversedRange; // bases covered by the cluster using ploidy data
    public long TotalDeleted; // bases with traversed range but without cluster CN support

    // output from chaining routine - % of copy number segments had discernable A, B and cluster ploidy
    public double ValidAllelePloidySegmentPerc;

    public ClusterMetrics()
    {
        DBCount = 0;
        ShortDBCount = 0;
        ClusterDBCount = 0;
        TotalDBLength = 0;
        TraversedDelLength = 0;
        TraversedDelCount = 0;
        ChainedLength = 0;
        TotalRange = 0;
        TraversedRange = 0;
        TotalDeleted = 0;
        ValidAllelePloidySegmentPerc = 1;
    }

}
