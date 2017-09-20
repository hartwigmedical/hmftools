package com.hartwig.hmftools.common.purple.ratio;

import java.util.List;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

class RollingRatioNormalization implements Supplier<List<ReadRatio>> {

    private int startIndex = 0;
    private int endIndex = -1;

    private final long maxWindowDistance;
    private final List<ReadRatio> ratios;
    private final List<ReadRatio> result = Lists.newArrayList();
    private final RollingMedian rollingMedian = new RollingMedian();

    RollingRatioNormalization(final double expectedRatio, final long maxWindowDistance, final long minWindowCoverage,
            final List<ReadRatio> ratios) {
        this.maxWindowDistance = maxWindowDistance;
        this.ratios = ratios;

        for (int currentIndex = 0; currentIndex < ratios.size(); currentIndex++) {
            final ReadRatio current = ratios.get(currentIndex);

            removeExpiredRatios(currentIndex);
            addNewRatios(currentIndex);

            double medianRatio = rollingMedian.median();
            double correctedRatio = isValid(current) && rollingMedian.size() >= minWindowCoverage
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

    private void addNewRatios(int currentIndex) {
        for (int laterIndex = endIndex + 1; laterIndex < ratios.size(); laterIndex++) {
            final ReadRatio later = ratios.get(laterIndex);

            if (distance(currentIndex, laterIndex) <= maxWindowDistance) {
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
        }
    }

    private void removeExpiredRatios(int currentIndex) {
        for (int earlierIndex = startIndex; earlierIndex < currentIndex; earlierIndex++) {
            final ReadRatio earlier = ratios.get(earlierIndex);
            final boolean isValid = isValid(earlier);

            if (!isValid || distance(currentIndex, earlierIndex) > maxWindowDistance) {
                if (isValid) {
                    rollingMedian.remove(earlier.ratio());
                }
                startIndex++;
            } else {
                return;
            }
        }
    }

    private long distance(int firstIndex, int secondIndex) {
        return Math.abs(firstIndex - secondIndex);
    }
}
