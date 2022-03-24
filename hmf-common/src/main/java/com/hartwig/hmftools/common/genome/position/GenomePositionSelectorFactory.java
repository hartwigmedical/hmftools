package com.hartwig.hmftools.common.genome.position;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public final class GenomePositionSelectorFactory {

    private GenomePositionSelectorFactory() {
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final Collection<P> positions) {
        return new GenomePositionSelectorIteratorImpl<>(positions);
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final List<P> positions) {
        return new GenomePositionSelectorListImpl<>(positions);
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final ListMultimap<Chromosome, P> positions) {
        final GenomePositionSelector<P> nullSelector = new NullGenomePositionSelector<>();

        final Map<Chromosome, GenomePositionSelector<P>> chromosomeSelectors = Maps.newHashMap();
        for (final Chromosome chromosome : positions.keySet()) {
            chromosomeSelectors.put(chromosome, new GenomePositionSelectorListImpl<>(positions.get(chromosome)));
        }
        return new GenomePositionSelector<P>() {
            @NotNull
            @Override
            public Optional<P> select(@NotNull final GenomePosition position) {
                return chromosomeSelectors.getOrDefault(HumanChromosome.fromString(position.chromosome()), nullSelector).select(position);
            }

            @Override
            public void select(final GenomeRegion region, final Consumer<P> handler) {
                chromosomeSelectors.getOrDefault(HumanChromosome.fromString(region.chromosome()), nullSelector).select(region, handler);
            }
        };
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final Map<Chromosome,List<P>> positions) {
        final GenomePositionSelector<P> nullSelector = new NullGenomePositionSelector<>();

        final Map<Chromosome, GenomePositionSelector<P>> chromosomeSelectors = Maps.newHashMap();
        for (final Chromosome chromosome : positions.keySet()) {
            chromosomeSelectors.put(chromosome, new GenomePositionSelectorListImpl<>(positions.get(chromosome)));
        }
        return new GenomePositionSelector<P>() {
            @NotNull
            @Override
            public Optional<P> select(@NotNull final GenomePosition position) {
                return chromosomeSelectors.getOrDefault(HumanChromosome.fromString(position.chromosome()), nullSelector).select(position);
            }

            @Override
            public void select(final GenomeRegion region, final Consumer<P> handler) {
                chromosomeSelectors.getOrDefault(HumanChromosome.fromString(region.chromosome()), nullSelector).select(region, handler);
            }
        };
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final Multimap<Chromosome, P> positions) {
        final GenomePositionSelector<P> nullSelector = new NullGenomePositionSelector<>();

        final Map<Chromosome, GenomePositionSelector<P>> chromosomeSelectors = Maps.newHashMap();
        for (final Chromosome chromosome : positions.keySet()) {
            chromosomeSelectors.put(chromosome, new GenomePositionSelectorIteratorImpl<>(positions.get(chromosome)));
        }
        return new GenomePositionSelector<P>() {
            @NotNull
            @Override
            public Optional<P> select(@NotNull final GenomePosition position) {
                return chromosomeSelectors.getOrDefault(HumanChromosome.fromString(position.chromosome()), nullSelector).select(position);
            }

            @Override
            public void select(final GenomeRegion region, final Consumer<P> handler) {
                chromosomeSelectors.getOrDefault(HumanChromosome.fromString(region.chromosome()), nullSelector).select(region, handler);
            }
        };
    }

    private static class NullGenomePositionSelector<P extends GenomePosition> implements GenomePositionSelector<P> {

        @NotNull
        @Override
        public Optional<P> select(@NotNull final GenomePosition position) {
            return Optional.empty();
        }

        @Override
        public void select(final GenomeRegion region, final Consumer<P> handler) {
            // VOID
        }
    }
}
