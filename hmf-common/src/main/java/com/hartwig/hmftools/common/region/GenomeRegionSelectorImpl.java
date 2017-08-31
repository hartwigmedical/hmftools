package com.hartwig.hmftools.common.region;

import java.util.Collection;
import java.util.Iterator;
import java.util.Optional;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomeRegionSelectorImpl<R extends GenomeRegion> implements GenomeRegionSelector<R> {

    @NotNull
    private final Iterator<R> regions;
    @Nullable
    private GenomePosition lastPosition;
    @Nullable
    private R next;

    GenomeRegionSelectorImpl(@NotNull Collection<R> regions) {
        this.regions = regions.iterator();
        next = this.regions.hasNext() ? this.regions.next() : null;
    }

    @NotNull
    public Optional<R> select(@NotNull final GenomePosition variant) {

        if (lastPosition != null && variant.compareTo(lastPosition) < 0) {
            throw new IllegalArgumentException("Selector only goes forward, never backwards!");
        }
        lastPosition = variant;
        while (next != null && compare(variant, next) > 0) {
            next = regions.hasNext() ? this.regions.next() : null;
        }

        return next != null && compare(variant, next) == 0 ? Optional.of(next) : Optional.empty();
    }

    public static int compare(@NotNull final GenomePosition position, @NotNull final GenomeRegion region) {
        int positionChromosome = HumanChromosome.fromString(position.chromosome()).intValue();
        int regionChromosome = HumanChromosome.fromString(region.chromosome()).intValue();
        if (positionChromosome < regionChromosome) {
            return -1;
        }
        if (positionChromosome > regionChromosome) {
            return 1;
        }

        if (position.position() < region.start()) {
            return -1;
        }

        if (position.position() > region.end()) {
            return 1;
        }

        return 0;
    }

}
