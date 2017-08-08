package com.hartwig.hmftools.common.position;

import static com.hartwig.hmftools.common.region.GenomeRegionSelectorImpl.compare;

import java.util.Collection;
import java.util.Iterator;
import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomePositionSelectorImpl<P extends GenomePosition> implements GenomePositionSelector<P> {

    @NotNull
    private final Iterator<P> positions;
    @Nullable
    private GenomeRegion lastRegion;
    @Nullable
    private P next;

    GenomePositionSelectorImpl(@NotNull Collection<P> positions) {
        this.positions = positions.iterator();
        next = this.positions.hasNext() ? this.positions.next() : null;
    }

    @Override
    public Optional<P> select(@NotNull final GenomePosition position) {

        while (next != null && next.compareTo(position) < 0) {
            next = positions.hasNext() ? this.positions.next() : null;
        }

        if (next != null && next.compareTo(position) == 0) {
            return Optional.of(next);
        }

        return Optional.empty();
    }

    @Override
    public void select(final GenomeRegion region, final Consumer<P> handler) {
        if (lastRegion != null && region.compareTo(lastRegion) < 0) {
            throw new IllegalArgumentException("Selector only goes forward, never backwards!");
        }
        lastRegion = region;

        while (next != null && compare(next, region) < 0) {
            next = positions.hasNext() ? this.positions.next() : null;
        }

        while (next != null && compare(next, region) == 0) {
            handler.accept(next);
            next = positions.hasNext() ? this.positions.next() : null;
        }
    }

}
