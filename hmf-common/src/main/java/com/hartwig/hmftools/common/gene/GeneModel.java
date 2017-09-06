package com.hartwig.hmftools.common.gene;

import java.util.Collection;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;

import org.jetbrains.annotations.NotNull;

public class GeneModel {

    private final Collection<HmfGenomeRegion> regions;
    private final Slicer slicer;

    public GeneModel(@NotNull final SortedSetMultimap<String, HmfGenomeRegion> regions) {
        this.regions = regions.values();
        slicer = SlicerFactory.fromRegions(regions);
    }

    public Collection<HmfGenomeRegion> hmfRegions() {
        return regions;
    }

    public Slicer slicer() {
        return slicer;
    }
}
