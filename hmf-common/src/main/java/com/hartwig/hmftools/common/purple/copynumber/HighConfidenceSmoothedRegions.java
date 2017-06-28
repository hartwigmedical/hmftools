package com.hartwig.hmftools.common.purple.copynumber;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class HighConfidenceSmoothedRegions extends BaseSmoothedRegions {

    private final double purity;
    private final List<PurpleCopyNumber> smoothedRegions = Lists.newArrayList();
    private final List<PurpleCopyNumber> highConfidenceRegions;
    private final List<FittedRegion> fittedRegions;

    HighConfidenceSmoothedRegions(double purity, @NotNull final List<PurpleCopyNumber> highConfidenceRegions,
            @NotNull final List<FittedRegion> fittedRegions) {
        this.purity = purity;
        this.highConfidenceRegions = highConfidenceRegions;
        this.fittedRegions = fittedRegions;

        run();
    }

    @NotNull
    List<PurpleCopyNumber> smoothedRegions() {
        return smoothedRegions;
    }

    private void run() {
        if (!highConfidenceRegions.isEmpty()) {
            int largestIncludedIndex = -1;
            HighConfidenceCopyNumberBuilder currentBuilder;

            for (int i = 0; i < highConfidenceRegions.size(); i++) {
                final PurpleCopyNumber currentRegion = highConfidenceRegions.get(i);
                int startOfRegionIndex = indexOfStart(largestIncludedIndex + 1, currentRegion);
                int endOfRegionIndex = indexOfEnd(startOfRegionIndex, currentRegion);

                // JOBA: Start new builder
                currentBuilder = new HighConfidenceCopyNumberBuilder(purity, fittedRegions.get(startOfRegionIndex));

                // JOBA: Go backwards to previous end
                currentBuilder = backwards(startOfRegionIndex - 1, largestIncludedIndex + 1, currentBuilder);

                // JOBA: Go forwards to end of region
                currentBuilder = forwards(startOfRegionIndex + 1, endOfRegionIndex, currentBuilder);

                boolean isLastBroadRegion = i == highConfidenceRegions.size() - 1;
                if (isLastBroadRegion) {
                    currentBuilder = forwards(endOfRegionIndex + 1, fittedRegions.size() - 1, currentBuilder);
                    smoothedRegions.add(currentBuilder.build());
                } else {
                    int nextStartIndex = indexOfStart(endOfRegionIndex + 1, highConfidenceRegions.get(i + 1));
                    largestIncludedIndex = forwardsUntilDifferent(endOfRegionIndex + 1, nextStartIndex - 1, currentBuilder);
                    smoothedRegions.add(currentBuilder.build());
                }
            }
        }
    }

    private int forwardsUntilDifferent(int startIndex, int endIndex, @NotNull HighConfidenceCopyNumberBuilder builder) {
        for (int i = startIndex; i <= endIndex; i++) {
            FittedRegion copyNumber = fittedRegions.get(i);
            if (isSimilar(copyNumber, builder)) {
                builder.extendRegion(copyNumber);
            } else {
                return i - 1;
            }
        }

        return endIndex;
    }

    @NotNull
    private HighConfidenceCopyNumberBuilder forwards(int startIndex, int endIndex, final @NotNull HighConfidenceCopyNumberBuilder builder) {
        HighConfidenceCopyNumberBuilder current = builder;
        for (int i = startIndex; i <= endIndex; i++) {
            FittedRegion copyNumber = fittedRegions.get(i);
            if (isSimilar(copyNumber, current)) {
                current.extendRegion(copyNumber);
            } else {
                smoothedRegions.add(current.build());
                current = new HighConfidenceCopyNumberBuilder(purity, copyNumber);
            }
        }

        return current;
    }

    @NotNull
    private HighConfidenceCopyNumberBuilder backwards(int startIndex, int endIndex,
            @NotNull final HighConfidenceCopyNumberBuilder forwardBuilder) {
        final Deque<PurpleCopyNumber> preRegions = new ArrayDeque<>();
        HighConfidenceCopyNumberBuilder reverseBuilder = forwardBuilder;

        for (int i = startIndex; i >= endIndex; i--) {
            final FittedRegion copyNumber = fittedRegions.get(i);
            if (isSimilar(copyNumber, reverseBuilder)) {
                reverseBuilder.extendRegion(copyNumber);
            } else {
                if (reverseBuilder != forwardBuilder) {
                    preRegions.addFirst(reverseBuilder.build());
                }
                reverseBuilder = new HighConfidenceCopyNumberBuilder(purity, copyNumber);
            }
        }

        if (reverseBuilder != forwardBuilder) {
            preRegions.addFirst(reverseBuilder.build());
        }

        smoothedRegions.addAll(preRegions);
        return forwardBuilder;
    }

    private static boolean isSimilar(@NotNull final FittedRegion copyNumber, @NotNull final HighConfidenceCopyNumberBuilder builder) {
        int bafCount = copyNumber.bafCount();
        if (!isDiploid(copyNumber)) {
            return true;
        }

        if (!builder.withinCopyNumberTolerance(copyNumber)) {
            return false;
        }

        if (bafCount > 0 && !Doubles.isZero(builder.averageObservedBAF())) {
            double bafDeviation = Math.abs(copyNumber.observedBAF() - builder.averageObservedBAF());
            if (Doubles.greaterThan(bafDeviation, allowedBAFDeviation(bafCount))) {
                return false;
            }
        }

        return true;
    }

    static double allowedBAFDeviation(int bafCount) {
        return 1d / Math.pow(Math.max(1, bafCount), 0.95) * 0.38 + 0.022;
    }

    private int indexOfEnd(int minIndex, @NotNull PurpleCopyNumber region) {
        return indexOf(minIndex, copyNumber -> copyNumber.end() == region.end());
    }

    private int indexOfStart(int minIndex, @NotNull PurpleCopyNumber region) {
        return indexOf(minIndex, copyNumber -> copyNumber.start() == region.start());
    }

    private int indexOf(int minIndex, @NotNull Predicate<FittedRegion> predicate) {
        for (int i = minIndex; i < fittedRegions.size(); i++) {
            FittedRegion copyNumber = fittedRegions.get(i);
            if (predicate.test(copyNumber)) {
                return i;
            }
        }

        throw new IllegalArgumentException();
    }
}
