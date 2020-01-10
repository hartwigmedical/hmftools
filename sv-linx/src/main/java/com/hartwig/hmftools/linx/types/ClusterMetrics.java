package com.hartwig.hmftools.linx.types;

public class ClusterMetrics
{
    public int DBCount;
    public int ShortDBCount;
    public int ClusterDBCount;
    public long TotalDeleted;
    public long ChainedLength;
    public long TotalRange;

    public ClusterMetrics()
    {
        DBCount = 0;
        ShortDBCount = 0;
        ClusterDBCount = 0;
        TotalDeleted = 0;
        ChainedLength = 0;
        TotalRange = 0;
    }

}
