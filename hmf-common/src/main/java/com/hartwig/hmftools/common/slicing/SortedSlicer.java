package com.hartwig.hmftools.common.slicing;

import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

import java.util.*;
import java.util.stream.Collectors;

class SortedSlicer implements Slicer {

    private final Map<String, SingleChromosomeSlicer> regions = Maps.newHashMap();

    SortedSlicer(@NotNull final SortedSetMultimap<String, GenomeRegion> regions) {
        for (final String chromosome : regions.keySet()) {
            SortedSet<GenomeRegion> region = regions.get(chromosome);
            this.regions.put(chromosome, new SingleChromosomeSlicer(region));
        }
    }

    @Override
    public boolean test(@NotNull final GenomePosition variant) {
        final SingleChromosomeSlicer slicer = regions.get(variant.chromosome());
        return slicer != null && slicer.includes(variant);
    }

    @NotNull
    @Override
    public Collection<GenomeRegion> regions() {
        return regions.values().stream().flatMap(x -> x.deque.stream()).collect(Collectors.toList());
    }

    private static class SingleChromosomeSlicer {

        @NotNull
        private final Deque<GenomeRegion> deque;
        private long currentPosition;

        SingleChromosomeSlicer(@NotNull SortedSet<GenomeRegion> region) {
            this.deque = new ArrayDeque<>(region);
        }

        boolean includes(@NotNull final GenomePosition variant) {
            if (variant.position() < currentPosition) {
                throw new IllegalArgumentException("Forward slicer only goes forward, never backwards!");
            }
            currentPosition = variant.position();

            GenomeRegion region = deque.peek();
            while (region != null && region.end() < variant.position()) {
                deque.pop();
                region = deque.peek();
            }
            return region != null && variant.position() >= region.start();
        }
    }
}
