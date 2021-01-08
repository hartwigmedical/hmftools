package com.hartwig.hmftools.common.genome.region;

import java.util.Collection;

import org.jetbrains.annotations.NotNull;

public class GenomeRegionsValidation {

    public static boolean isSubset(@NotNull Collection<? extends GenomeRegion> superset, @NotNull Collection<? extends GenomeRegion> subset) {
        GenomeRegionSelector<? extends GenomeRegion> selector = GenomeRegionSelectorFactory.create(superset);
        for (GenomeRegion region : subset) {
            GenomeRegionsBuilder builder = new GenomeRegionsBuilder();
            selector.select(region, builder::addRegion);
            if (builder.build().stream().noneMatch(x -> x.overlappingBases(region) == region.bases())) {
                return false;
            }
        }

        return true;
    }

}
