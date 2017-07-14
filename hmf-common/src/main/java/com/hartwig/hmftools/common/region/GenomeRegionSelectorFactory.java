package com.hartwig.hmftools.common.region;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public class GenomeRegionSelectorFactory {

    @NotNull
    public static <R extends GenomeRegion> GenomeRegionSelector<R> create(@NotNull final List<R> regions) {
        return new GenomeRegionSelectorImpl<>(regions);
    }

    @NotNull
    public static <R extends GenomeRegion> GenomeRegionSelector<R> create(@NotNull final Multimap<String, R> regions) {

        final GenomeRegionSelector<R> nullSelector = new NullGenomeRegionSelector<>();

        final Map<String, GenomeRegionSelector<R>> chromosomeSelectors = Maps.newHashMap();
        for (final String chromosome : regions.keySet()) {
            chromosomeSelectors.put(chromosome, new GenomeRegionSelectorImpl<>(regions.get(chromosome)));
        }

        return position -> chromosomeSelectors.getOrDefault(position.chromosome(), nullSelector).select(position);
    }

    private static class NullGenomeRegionSelector<R extends GenomeRegion> implements GenomeRegionSelector<R> {

        @NotNull
        @Override
        public Optional<R> select(final GenomePosition position) {
            return Optional.empty();
        }
    }
}
