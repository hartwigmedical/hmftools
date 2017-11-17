package com.hartwig.hmftools.common.region;

import java.util.Collection;
import java.util.Iterator;
import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

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

    @Override
    public void select(@NotNull final GenomeRegion region, @NotNull final Consumer<R> handler) {
        final GenomePosition start = GenomePositions.create(region.chromosome(), region.start());
        final GenomePosition end = GenomePositions.create(region.chromosome(), region.end());
        if (lastPosition != null && start.compareTo(lastPosition) < 0) {
            throw new IllegalArgumentException("Selector only goes forward, never backwards!");
        }

        while (next != null && compare(start, next) > 0) {
            next = regions.hasNext() ? this.regions.next() : null;
        }

        while (next != null && compare(start, next) <= 0 && compare(end, next) >= 0 ) {
            handler.accept(next);
            next = regions.hasNext() ? this.regions.next() : null;
        }

        lastPosition = end;
    }

    @NotNull
    public Optional<R> select(@NotNull final GenomePosition position) {

        if (lastPosition != null && position.compareTo(lastPosition) < 0) {
            throw new IllegalArgumentException("Selector only goes forward, never backwards!");
        }
        lastPosition = position;
        while (next != null && compare(position, next) > 0) {
            next = regions.hasNext() ? this.regions.next() : null;
        }

        return next != null && compare(position, next) == 0 ? Optional.of(next) : Optional.empty();
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
