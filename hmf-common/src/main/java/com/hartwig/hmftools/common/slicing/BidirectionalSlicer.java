package com.hartwig.hmftools.common.slicing;

import java.util.Collection;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

class BidirectionalSlicer implements Slicer {

    @NotNull
    private final Multimap<String, GenomeRegion> regions;

    BidirectionalSlicer(@NotNull final Multimap<String, GenomeRegion> regions) {
        this.regions = regions;
    }

    @Override
    public boolean test(@NotNull GenomePosition variant) {
        final Collection<GenomeRegion> regionsForChrom = regions.get(variant.chromosome());
        if (regionsForChrom == null) {
            return false;
        } else {
            for (final GenomeRegion region : regionsForChrom) {
                if (variant.position() >= region.start() && variant.position() <= region.end()) {
                    return true;
                } else if (region.start() > variant.position()) {
                    return false;
                }
            }
        }

        return false;
    }

    @NotNull
    @Override
    public Collection<GenomeRegion> regions() {
        return regions.values();
    }
}
