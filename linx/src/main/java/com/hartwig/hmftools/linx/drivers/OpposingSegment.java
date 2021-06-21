package com.hartwig.hmftools.linx.drivers;

import static java.lang.Math.max;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;

public class OpposingSegment
{
    public final SvCluster Cluster;
    public final List<SvBreakend> Breakends;
    public final double NetCNChange;

    private double mRemainingCNChange;

    public OpposingSegment(final SvCluster cluster, final List<SvBreakend> breakends, double netCNChange)
    {
        Cluster = cluster;
        Breakends = breakends;
        NetCNChange = netCNChange;
        mRemainingCNChange = netCNChange;
    }

    public double remainingCNChange() { return mRemainingCNChange; }
    public void zeroCNChange() { mRemainingCNChange = 0; }
    public void reduceCNChange(double change) { mRemainingCNChange = max(mRemainingCNChange - change, 0); }

    public String toString()
    {
        return String.format("cluster(%d) breakends(%d) cnChange(orig=%.1f remaining=%.1f)",
                Cluster.id(), Breakends.size(), NetCNChange, mRemainingCNChange);
    }

}
