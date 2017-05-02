package com.hartwig.hmftools.common.slicing;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.variant.Variant;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;

class UnsortedSlicer implements Slicer {

    @NotNull
    private final SortedSetMultimap<String, GenomeRegion> regions;

    UnsortedSlicer(@NotNull final SortedSetMultimap<String, GenomeRegion> regions) {
        this.regions = regions;
    }

    public boolean test(@NotNull Variant variant) {
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
