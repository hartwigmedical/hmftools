package com.hartwig.hmftools.common.region;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomeRegionSelector<R extends GenomeRegion> {

    @NotNull
    private final Map<String, SingleChromosomeSelector<R>> chromosomeSelectors = Maps.newHashMap();

    public GenomeRegionSelector(@NotNull final List<R> regions) {
        this(Multimaps.index(regions, GenomeRegion::chromosome));
    }

    public GenomeRegionSelector(@NotNull final Multimap<String, R> regions) {
        for (final String chromosome : regions.keySet()) {
            SingleChromosomeSelector<R> selector = new SingleChromosomeSelector<>(regions.get(chromosome));
            this.chromosomeSelectors.put(chromosome, selector);
        }
    }

    @NotNull
    public Optional<R> select(GenomePosition position) {
        return Optional.ofNullable(chromosomeSelectors.get(position.chromosome())).flatMap(x -> x.region(position));
    }

    private static class SingleChromosomeSelector<R extends GenomeRegion> {

        @NotNull
        private final Iterator<R> regions;
        private long currentPosition;
        @Nullable
        private R next;

        SingleChromosomeSelector(@NotNull Collection<R> regions) {
            this.regions = regions.iterator();
            next = this.regions.hasNext() ? this.regions.next() : null;
        }

        Optional<R> region(@NotNull final GenomePosition variant) {
            if (variant.position() < currentPosition) {
                throw new IllegalArgumentException("Selector only goes forward, never backwards!");
            }
            currentPosition = variant.position();
            while (next != null && next.end() < variant.position()) {
                next = regions.hasNext() ? this.regions.next() : null;
            }

            return next != null && variant.position() >= next.start() ? Optional.of(next) : Optional.empty();
        }
    }
}
