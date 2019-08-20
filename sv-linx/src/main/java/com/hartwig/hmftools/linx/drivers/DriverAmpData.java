package com.hartwig.hmftools.linx.drivers;

import com.hartwig.hmftools.linx.types.SvCluster;

public class DriverAmpData
{
    public final SvCluster Cluster;

    public final boolean TraverseUp;
    public final int BreakendCount;
    public final int SegmentCount;
    public final double StartCopyNumber;
    public final double NetCNChange;

    public DriverAmpData(final SvCluster cluster, boolean traverseUp, int breakendCount, int segmentCount,
            double startCopyNumber, double cnChange)
    {
        Cluster = cluster;
        TraverseUp = traverseUp;
        BreakendCount = breakendCount;
        SegmentCount = segmentCount;
        StartCopyNumber = startCopyNumber;
        NetCNChange = cnChange;
    }

    public String toString()
    {
        return String.format("cluster(%d) netCNChange(%.1f) startCN(%.1f) segs(%d) breakends(%d) traversal(%s)",
                Cluster.id(), NetCNChange, StartCopyNumber, SegmentCount, BreakendCount, TraverseUp ? "up" : "down");
    }

}
