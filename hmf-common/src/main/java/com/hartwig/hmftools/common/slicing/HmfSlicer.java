package com.hartwig.hmftools.common.slicing;

import java.util.Collection;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.jetbrains.annotations.NotNull;

public class HmfSlicer extends BidirectionalSlicer {

    @NotNull
    private final SortedSetMultimap<String, HmfGenomeRegion> regions;

    HmfSlicer(@NotNull final SortedSetMultimap<String, HmfGenomeRegion> regions) {
        super(TreeMultimap.create(regions));
        this.regions = regions;
    }

    @NotNull
    public Collection<HmfGenomeRegion> hmfRegions() {
        return regions.values();
    }
}
