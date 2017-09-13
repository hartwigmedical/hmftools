package com.hartwig.hmftools.common.purple.ratio;

import java.util.List;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

class RollingRatioNormalization implements Supplier<List<ReadRatio>> {

    private int startIndex = 0;
    private int endIndex = -1;

    private final long distance;
    private final List<ReadRatio> result = Lists.newArrayList();
    private final RollingMedian rollingMedian = new RollingMedian();

    RollingRatioNormalization(final double expectedRatio, final long distance, final List<ReadRatio> ratios) {
        this.distance = distance;

        initializeRollingMedian(ratios);

        for (final ReadRatio current : ratios) {
            removeExpiredRatios(current, ratios);
            addNewRatios(current, ratios);

            double medianRatio = rollingMedian.medianRatio();
            double correctedRatio = isValid(current) ? expectedRatio * current.ratio() / medianRatio : -1;

            result.add(ImmutableReadRatio.builder().from(current).ratio(correctedRatio).build());
        }
    }

    private boolean isValid(@NotNull final ReadRatio ratio) {
        return !Doubles.equal(ratio.ratio(), -1);
    }

    @Override
    public List<ReadRatio> get() {
        return result;
    }

    private void initializeRollingMedian(@NotNull final List<ReadRatio> ratios) {
        startIndex = 0;
        for (ReadRatio ratio : ratios) {
            if (ratio.position() <= distance) {
                addToMedian(ratio);
            } else {
                return;
            }
        }
    }

    private void removeExpiredRatios(@NotNull final ReadRatio current, @NotNull final List<ReadRatio> ratios) {
        for (int earlierIndex = startIndex; earlierIndex < ratios.size(); earlierIndex++) {
            final ReadRatio earlier = ratios.get(earlierIndex);

            if (distance(current, earlier) > distance) {
                if (isValid(earlier)) {
                    rollingMedian.remove(earlier.ratio());
                }
                startIndex++;
            } else {
                return;
            }
        }
    }

    private void addNewRatios(@NotNull final ReadRatio current, @NotNull final List<ReadRatio> ratios) {
        for (int i = endIndex + 1; i < ratios.size(); i++) {
            final ReadRatio later = ratios.get(i);

            if (distance(current, later) <= distance) {
                addToMedian(later);
            } else {
                return;
            }
        }
    }

    private void addToMedian(@NotNull final ReadRatio current) {
        if (isValid(current)) {
            rollingMedian.add(current.ratio());
        }
        endIndex++;
    }

    private long distance(@NotNull final ReadRatio first, @NotNull final ReadRatio second) {
        return Math.abs(first.position() - second.position());
    }
}
