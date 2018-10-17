package com.hartwig.hmftools.common.region;

import java.util.Collection;
import java.util.Map;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public final class GenomeRegionSelectorFactory {

    private GenomeRegionSelectorFactory() {
    }

    @NotNull
    public static <R extends GenomeRegion> GenomeRegionSelector<R> create(@NotNull final Collection<R> regions) {
        return new GenomeRegionSelectorImpl<>(regions);
    }

    @NotNull
    public static <R extends GenomeRegion> GenomeRegionSelector<R> create(@NotNull final Multimap<String, R> regions) {
        final GenomeRegionSelector<R> nullSelector = new NullGenomeRegionSelector<>();

        final Map<String, GenomeRegionSelector<R>> chromosomeSelectors = Maps.newHashMap();
        for (final String chromosome : regions.keySet()) {
            chromosomeSelectors.put(chromosome, new GenomeRegionSelectorImpl<>(regions.get(chromosome)));
        }

        return new GenomeRegionSelector<R>() {
            @NotNull
            @Override
            public Optional<R> select(@NotNull final GenomePosition position) {
                return chromosomeSelectors.getOrDefault(position.chromosome(), nullSelector).select(position);
            }

            @Override
            public void select(@NotNull final GenomeRegion region, @NotNull final Consumer<R> handler) {
                chromosomeSelectors.getOrDefault(region.chromosome(), nullSelector).select(region, handler);
            }
        };
    }

    @NotNull
    public static <R extends GenomeRegion> GenomeRegionSelector<R> createImproved(@NotNull final Multimap<Chromosome, R> regions) {
        final GenomeRegionSelector<R> nullSelector = new NullGenomeRegionSelector<>();

        final Map<Chromosome, GenomeRegionSelector<R>> chromosomeSelectors = Maps.newHashMap();
        for (final Chromosome chromosome : regions.keySet()) {
            chromosomeSelectors.put(chromosome, new GenomeRegionSelectorImpl<>(regions.get(chromosome)));
        }

        return new GenomeRegionSelector<R>() {
            @NotNull
            @Override
            public Optional<R> select(@NotNull final GenomePosition position) {
                return chromosomeSelectors.getOrDefault(HumanChromosome.fromString(position.chromosome()), nullSelector).select(position);
            }

            @Override
            public void select(@NotNull final GenomeRegion region, @NotNull final Consumer<R> handler) {
                chromosomeSelectors.getOrDefault(HumanChromosome.fromString(region.chromosome()), nullSelector).select(region, handler);
            }
        };
    }

    private static class NullGenomeRegionSelector<R extends GenomeRegion> implements GenomeRegionSelector<R> {

        @NotNull
        @Override
        public Optional<R> select(@NotNull final GenomePosition position) {
            return Optional.empty();
        }

        @Override
        public void select(@NotNull final GenomeRegion region, @NotNull final Consumer<R> handler) {

        }
    }
}
