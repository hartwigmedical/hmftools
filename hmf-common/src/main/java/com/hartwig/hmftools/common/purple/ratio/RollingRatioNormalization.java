package com.hartwig.hmftools.common.purple.ratio;

import java.util.List;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

class RollingRatioNormalization implements Supplier<List<ReadRatio>> {

    private int startIndex = 0;
    private int endIndex = -1;
    private int endValidIndex = -1;

    private final long maxDistance;
    private final long minCoverage;
    private final List<ReadRatio> ratios;
    private final List<ReadRatio> result = Lists.newArrayList();
    private final RollingMedian rollingMedian = new RollingMedian();

    RollingRatioNormalization(final double expectedRatio, final long maxDistance, final long minCoverage, final List<ReadRatio> ratios) {
        this.maxDistance = maxDistance;
        this.minCoverage = minCoverage;
        this.ratios = ratios;

        initializeRollingMedian(ratios);

        for (int currentIndex = 0; currentIndex < ratios.size(); currentIndex++) {
            final ReadRatio current = ratios.get(currentIndex);

            removeExpiredRatios(currentIndex);
            addNewRatios(currentIndex);

            double medianRatio = rollingMedian.median();
            double correctedRatio = isValid(current) && sufficientCoverage() && !rollingMedian.isEmpty()
                    ? expectedRatio * current.ratio() / medianRatio
                    : current.ratio();

            result.add(ImmutableReadRatio.builder().from(current).ratio(correctedRatio).build());
        }
    }

    @Override
    public List<ReadRatio> get() {
        return result;
    }

    private boolean isValid(@NotNull final ReadRatio ratio) {
        return Doubles.greaterThan(ratio.ratio(), 0);
    }


    private void initializeRollingMedian(@NotNull final List<ReadRatio> ratios) {
        startIndex = 0;
        for (ReadRatio ratio : ratios) {
            if (ratio.position() <= maxDistance) {
                addToMedian(ratio);
            } else {
                return;
            }
        }
    }

    private void addNewRatios(int currentIndex) {
        final ReadRatio current = ratios.get(currentIndex);
        for (int i = endIndex + 1; i < ratios.size(); i++) {
            final ReadRatio later = ratios.get(i);

            if (distance(current, later) <= maxDistance) {
                addToMedian(later);
            } else {
                return;
            }
        }
    }

    private void addToMedian(@NotNull final ReadRatio current) {
        endIndex++;
        if (isValid(current)) {
            rollingMedian.add(current.ratio());
            endValidIndex = endIndex;
        }
    }

    private void removeExpiredRatios(int currentIndex) {
        final ReadRatio current = ratios.get(currentIndex);
        for (int earlierIndex = startIndex; earlierIndex < currentIndex; earlierIndex++) {
            final ReadRatio earlier = ratios.get(earlierIndex);
            final boolean isValid = isValid(earlier);

            if (!isValid || distance(current, earlier) > maxDistance) {
                if (isValid) {
                    rollingMedian.remove(earlier.ratio());
                }
                startIndex++;
            } else {
                return;
            }
        }
    }

    private boolean sufficientCoverage() {
        return endValidIndex > startIndex && ratios.get(endValidIndex).position() - ratios.get(startIndex).position() >= minCoverage;
    }

    private long distance(@NotNull final ReadRatio first, @NotNull final ReadRatio second) {
        return Math.abs(first.position() - second.position());
    }
}
