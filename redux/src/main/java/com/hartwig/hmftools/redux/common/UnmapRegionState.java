package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class UnmapRegionState
{
    public final ChrBaseRegion Partition;
    public final List<HighDepthRegion> PartitionRegions;
    public Integer LastMatchedRegionIndex;

    public UnmapRegionState(final ChrBaseRegion partition, final List<HighDepthRegion> partitionRegions)
    {
        Partition = partition;
        PartitionRegions = partitionRegions;
        LastMatchedRegionIndex = null;
    }

    public String toString() { return format("partition(%s) regions(%d) index(%d)",
            Partition, PartitionRegions.size(), LastMatchedRegionIndex != null ? LastMatchedRegionIndex : -1); }
}
