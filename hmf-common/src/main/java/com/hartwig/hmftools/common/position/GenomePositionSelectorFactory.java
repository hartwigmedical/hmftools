package com.hartwig.hmftools.common.position;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class GenomePositionSelectorFactory {

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final List<P> positions) {
        return new GenomePositionSelectorImpl<>(positions);
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final Multimap<String, P> positions) {

        final GenomePositionSelector<P> nullSelector = new NullGenomePositionSelector<>();

        final Map<String, GenomePositionSelector<P>> chromosomeSelectors = Maps.newHashMap();
        for (final String chromosome : positions.keySet()) {
            chromosomeSelectors.put(chromosome, new GenomePositionSelectorImpl<>(positions.get(chromosome)));
        }
        return new GenomePositionSelector<P>() {
            @NotNull
            @Override
            public Optional<P> select(@NotNull final GenomePosition position) {
                return chromosomeSelectors.getOrDefault(position.chromosome(), nullSelector).select(position);
            }

            @Override
            public void select(final GenomeRegion region, final Consumer<P> handler) {
                chromosomeSelectors.getOrDefault(region.chromosome(), nullSelector).select(region, handler);
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
