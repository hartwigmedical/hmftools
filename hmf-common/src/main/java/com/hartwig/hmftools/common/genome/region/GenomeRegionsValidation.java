package com.hartwig.hmftools.common.genome.region;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class GenomeRegionsValidation {

    public static boolean isSubset(@NotNull Collection<? extends GenomeRegion> superset,
            @NotNull Collection<? extends GenomeRegion> subset) {
        final List<GenomeRegion> supersetList = superset.stream().sorted().collect(Collectors.toList());
        final List<GenomeRegion> subsetList = subset.stream().sorted().collect(Collectors.toList());

        GenomeRegionSelector<? extends GenomeRegion> selector = GenomeRegionSelectorFactory.create(supersetList);
        for (GenomeRegion region : subsetList) {
            GenomeRegionsBuilder builder = new GenomeRegionsBuilder();
            selector.select(region, builder::addRegion);
            if (builder.build().stream().noneMatch(x -> x.overlappingBases(region) == region.bases())) {
                return false;
            }
        }

        return true;
    }

}
