package com.hartwig.hmftools.common.slicing;


import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.variant.Variant;
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

    public boolean test(@NotNull Variant variant) {
        SingleChromosomeSlicer slicer = regions.get(variant.chromosome());
        return slicer != null && slicer.includes(variant);
    }

    @NotNull
    @Override
    public Collection<GenomeRegion> regions() {
        return regions.values().stream().flatMap(x -> x.deque.stream()).collect(Collectors.toList());
    }

    class SingleChromosomeSlicer {

        private final Deque<GenomeRegion> deque;
        private long currentPosition;

        SingleChromosomeSlicer(SortedSet<GenomeRegion> reqion) {
            this.deque = new ArrayDeque<>(reqion);
        }

        boolean includes(@NotNull Variant variant) {
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
