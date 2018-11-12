package com.hartwig.hmftools.common.position;

import static com.hartwig.hmftools.common.region.GenomeRegionSelectorImpl.compare;

import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class GenomePositionSelectorListImpl<P extends GenomePosition> implements GenomePositionSelector<P> {

    @NotNull
    private final List<P> positions;
    private int index = 0;

    GenomePositionSelectorListImpl(@NotNull List<P> positions) {
        this.positions = positions;
    }

    @Override
    @NotNull
    public Optional<P> select(@NotNull final GenomePosition position) {
        if (positions.isEmpty()) {
            return Optional.empty();
        }

        int currentCompare = current().compareTo(position);
        while (currentCompare > 0 && index > 0) {
            index--;
            currentCompare = current().compareTo(position);
        }

        while (currentCompare < 0 && index < positions.size() - 1) {
            index++;
            currentCompare = current().compareTo(position);
        }

        return currentCompare == 0 ? Optional.of(current()) : Optional.empty();
    }

    @Override
    public void select(final GenomeRegion region, final Consumer<P> handler) {
        if (positions.isEmpty()) {
            return;
        }

        // Line up index to start
        select(GenomePositions.create(region.chromosome(), region.start()));

        while (index < positions.size() && compare(current(), region) == 0) {
            handler.accept(current());
            index++;
        }

        index = Math.min(index, positions.size() - 1);
    }

    @NotNull
    private P current() {
        return positions.get(index);
    }

}
