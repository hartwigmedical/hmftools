package com.hartwig.hmftools.common.slicing;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ForwardSlicer implements Slicer {

    private final Multimap<String, GenomeRegion> allRegions;
    private final Map<String, SingleChromosomeSlicer> chromosomeSlicers = Maps.newHashMap();

    ForwardSlicer(@NotNull final Multimap<String, GenomeRegion> regions) {
        allRegions = regions;
        for (final String chromosome : regions.keySet()) {
            this.chromosomeSlicers.put(chromosome, new SingleChromosomeSlicer(regions.get(chromosome)));
        }
    }

    @Override
    public boolean test(@NotNull final GenomePosition variant) {
        final SingleChromosomeSlicer slicer = chromosomeSlicers.get(variant.chromosome());
        return slicer != null && slicer.includes(variant);
    }

    @NotNull
    @Override
    public Collection<GenomeRegion> regions() {
        return allRegions.values();
    }

    private static class SingleChromosomeSlicer {

        @NotNull
        private final Iterator<GenomeRegion> regions;
        private long currentPosition;
        @Nullable
        private GenomeRegion next;

        SingleChromosomeSlicer(@NotNull Collection<GenomeRegion> regions) {
            this.regions = regions.iterator();
            next = this.regions.hasNext() ? this.regions.next() : null;
        }

        boolean includes(@NotNull final GenomePosition variant) {
            if (variant.position() < currentPosition) {
                throw new IllegalArgumentException("Forward slicer only goes forward, never backwards!");
            }
            currentPosition = variant.position();
            while (next != null && next.end() < variant.position()) {
                next = regions.hasNext() ? this.regions.next() : null;
            }

            return next != null && variant.position() >= next.start();
        }
    }
}
