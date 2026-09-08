package com.hartwig.hmftools.bamtools.depth;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Integers;

class PositionCount
{
    public int Position;
    public int DepthMin;
    public int DepthMax;
    public final List<Integer> DepthAverages;

    public PositionCount(final int position, int depthMin, int depthMax, int depthAvg)
    {
        Position = position;
        DepthMax = depthMax;
        DepthMin = depthMin;
        DepthAverages = Lists.newArrayList(depthAvg);
    }

    public int count() { return DepthAverages.size(); }
    public int medianDepth() { return (int)Math.round(Integers.median(DepthAverages)); }

    public String toString() { return format("%d: depth range(%d-%d) samples(%d)", Position, DepthMin, DepthMax, DepthAverages.size()); }
}
