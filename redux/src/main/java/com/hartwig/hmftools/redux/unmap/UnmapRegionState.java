package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.UnmappingRegion;

public class UnmapRegionState
{
    public final ChrBaseRegion SliceRegion;
    public final List<UnmappingRegion> Regions;
    public Integer LastMatchedRegionIndex;

    public UnmapRegionState(final ChrBaseRegion sliceRegion, final List<UnmappingRegion> regions)
    {
        SliceRegion = sliceRegion;
        Regions = regions;
        LastMatchedRegionIndex = null;
    }

    public String toString() { return format("partition(%s) regions(%d) index(%d)",
            SliceRegion, Regions.size(), LastMatchedRegionIndex != null ? LastMatchedRegionIndex : -1); }
}
