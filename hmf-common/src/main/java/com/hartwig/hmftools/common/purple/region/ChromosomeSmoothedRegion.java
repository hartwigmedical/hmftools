package com.hartwig.hmftools.common.purple.region;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.FittedCopyNumber;

public class ChromosomeSmoothedRegion {

    private final List<ConsolidatedRegion> smoothedRegions = Lists.newArrayList();

    private final String chromosome;
    private final List<ConsolidatedRegion> broadRegions;
    private final List<FittedCopyNumber> copyNumbers;

    private static final double DIPLOID_MIN_RATIO = 0.75;
    private static final double DIPLOID_MAX_RATIO = 1.25;

    private static final int MAX_BAF_COUNT = 50;
    private static final double MAX_BAF_RANGE = 0.12;
    private static final double MIN_BAF_RANGE = 0.03;

    private static final double MIN_RATIO_RANGE = 0.25;
    private static final double MAX_RATIO_RANGE = 1.25;

    public ChromosomeSmoothedRegion(final String chromosome, final List<ConsolidatedRegion> broadRegions,
            final List<FittedCopyNumber> copyNumbers) {
        this.chromosome = chromosome;
        this.broadRegions = broadRegions;
        this.copyNumbers = copyNumbers;

        doStuff();
    }

    public List<ConsolidatedRegion> getSmoothedRegions() {
        return smoothedRegions;
    }

    private void doStuff() {

        if (!broadRegions.isEmpty()) {

            int largestIncludedIndex = -1;
            ConsolidatedRegionBuilder currentBuilder = null;

            for (int i = 0; i < broadRegions.size(); i++) {

                ConsolidatedRegion currentRegion = broadRegions.get(i);
                int startOfRegionIndex = indexOfStart(largestIncludedIndex + 1, currentRegion);
                int endOfRegionIndex = indexOfEnd(startOfRegionIndex, currentRegion);

                // Start new builder
                currentBuilder = new ConsolidatedRegionBuilder(copyNumbers.get(startOfRegionIndex));

                // Go backwards to previous end (+1)
                currentBuilder = backwards(startOfRegionIndex - 1, largestIncludedIndex + 1, currentBuilder);

                // Go forwards to end of region
                currentBuilder = forwards(startOfRegionIndex + 1, endOfRegionIndex, currentBuilder);

                boolean isLastBroadRegion = i == broadRegions.size() - 1;
                if (isLastBroadRegion) {
                    currentBuilder = forwards(endOfRegionIndex + 1, copyNumbers.size() - 1, currentBuilder);
                    smoothedRegions.add(currentBuilder.build());
                } else {
                    int nextStartIndex = indexOfStart(endOfRegionIndex + 1, broadRegions.get(i + 1));
                    largestIncludedIndex = forwardsUntilDifferent(endOfRegionIndex + 1, nextStartIndex - 1, currentBuilder);
                    smoothedRegions.add(currentBuilder.build());
                }

            }
        }
    }

    private int forwardsUntilDifferent(int startIndex, int endIndex, ConsolidatedRegionBuilder builder) {
        for (int i = startIndex; i <= endIndex; i++) {
            FittedCopyNumber copyNumber = copyNumbers.get(i);
            if (isAlmostSimilar(copyNumber, builder)) {
                builder.extendRegion(copyNumber);
            } else {
                return i - 1;
            }
        }

        return endIndex;
    }

    private ConsolidatedRegionBuilder forwards(int startIndex, int endIndex, ConsolidatedRegionBuilder builder) {
        for (int i = startIndex; i <= endIndex; i++) {
            FittedCopyNumber copyNumber = copyNumbers.get(i);
            if (isVerySimilar(copyNumber, builder)) {
                builder.extendRegion(copyNumber);
            } else {
                smoothedRegions.add(builder.build());
                builder = new ConsolidatedRegionBuilder(copyNumber);
            }
        }

        return builder;
    }

    private ConsolidatedRegionBuilder backwards(int startIndex, int endIndex, ConsolidatedRegionBuilder forwardBuilder) {

        final Deque<ConsolidatedRegion> preRegions = new ArrayDeque<>();
        ConsolidatedRegionBuilder reverseBuilder = forwardBuilder;

        for (int i = startIndex; i >= endIndex; i--) {

            FittedCopyNumber copyNumber = copyNumbers.get(i);
            if (isVerySimilar(copyNumber, reverseBuilder)) {
                reverseBuilder.extendRegion(copyNumber);
            } else {
                if (reverseBuilder != forwardBuilder) {
                    preRegions.addFirst(reverseBuilder.build());
                }
                reverseBuilder = new ConsolidatedRegionBuilder(copyNumber);
            }
        }

        // Finalise
        if (reverseBuilder != forwardBuilder) {
            preRegions.addFirst(reverseBuilder.build());
        }

        smoothedRegions.addAll(preRegions);
        return forwardBuilder;
    }

    private boolean isVerySimilar(FittedCopyNumber copyNumber, ConsolidatedRegionBuilder builder) {
        return isSimilar(copyNumber, builder);
    }

    private boolean isAlmostSimilar(FittedCopyNumber copyNumber, ConsolidatedRegionBuilder builder) {
        return isSimilar(copyNumber, builder);
    }

    private boolean isSimilar(FittedCopyNumber copyNumber, ConsolidatedRegionBuilder builder) {

        int bafCount = copyNumber.bafCount();
        if (!isDiploid(copyNumber)) {
            return true;
        }

        double ratioOfRatioDeviation = Math.abs(copyNumber.ratioOfRatios() - builder.averageRatioOfRatios());
        double tumorRatioDeviation = Math.abs(copyNumber.normalisedTumorRatio() - builder.averageRatioOfRatios());
        double ratioDeviation = Math.min(ratioOfRatioDeviation, tumorRatioDeviation);
        if (Doubles.greaterThan(ratioDeviation, allowedRatioDeviation(bafCount))) {
            return false;
        }

        if (bafCount > 0 && !Doubles.isZero(builder.averageBAF())) {
            double bafDeviation = Math.abs(copyNumber.actualBAF() - builder.averageBAF());
            if (Doubles.greaterThan(bafDeviation, allowedBAFDeviation(bafCount))) {
                return false;
            }
        }


        return true;
    }

    static double allowedBAFDeviation(int bafCount) {
        if (bafCount >= MAX_BAF_COUNT) {
            return MIN_BAF_RANGE;
        }

        return (MIN_BAF_RANGE - MAX_BAF_RANGE) / MAX_BAF_COUNT * bafCount + MAX_BAF_RANGE;
    }

    static double allowedRatioDeviation(int bafCount) {
        if (bafCount >= 10) {
            return MIN_RATIO_RANGE;
        }

        return (MIN_RATIO_RANGE - MAX_RATIO_RANGE) / 10 * bafCount + MAX_RATIO_RANGE;
    }

    private boolean isDiploid(FittedCopyNumber copyNumber) {
        return Doubles.greaterOrEqual(copyNumber.normalCNVRatio(), DIPLOID_MIN_RATIO) && Doubles.lessOrEqual(
                copyNumber.normalCNVRatio(), DIPLOID_MAX_RATIO);
    }

    private int indexOfEnd(int minIndex, ConsolidatedRegion region) {
        return indexOf(minIndex, copyNumber -> copyNumber.end() == region.end());
    }

    private int indexOfStart(int minIndex, ConsolidatedRegion region) {
        return indexOf(minIndex, copyNumber -> copyNumber.start() == region.start());
    }

    private int indexOf(int minIndex, Predicate<FittedCopyNumber> predicate) {
        for (int i = minIndex; i < copyNumbers.size(); i++) {
            FittedCopyNumber copyNumber = copyNumbers.get(i);
            if (predicate.test(copyNumber)) {
                return i;
            }
        }

        throw new IllegalArgumentException();
    }

}
