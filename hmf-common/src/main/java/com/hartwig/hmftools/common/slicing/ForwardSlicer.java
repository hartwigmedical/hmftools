package com.hartwig.hmftools.common.slicing;

import java.util.Collection;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;

import org.jetbrains.annotations.NotNull;

class ForwardSlicer implements Slicer {

    @NotNull
    final Multimap<String, GenomeRegion> regions;

    @NotNull
    private final GenomeRegionSelector<GenomeRegion> selector;

    ForwardSlicer(@NotNull final Multimap<String, GenomeRegion> regions) {
        this.regions = regions;
        this.selector = new GenomeRegionSelector<>(regions);
    }

    @NotNull
    @Override
    public Collection<GenomeRegion> regions() {
        return regions.values();
    }

    @Override
    public boolean test(final GenomePosition genomePosition) {
        return selector.select(genomePosition).isPresent();
    }
}
