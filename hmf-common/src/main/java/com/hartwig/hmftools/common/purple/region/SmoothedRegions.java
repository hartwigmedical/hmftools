package com.hartwig.hmftools.common.purple.region;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.FittedCopyNumber;

import org.jetbrains.annotations.Nullable;

class SmoothedRegions {

    private static final double DIPLOID_MIN_RATIO = 0.75;
    private static final double DIPLOID_MAX_RATIO = 1.25;

    private static final int MAX_BAF_COUNT = 50;
    private static final double MAX_BAF_RANGE = 0.12;
    private static final double MIN_BAF_RANGE = 0.03;
    private static final double RATIO_RANGE = 0.25;

    private final List<ConsolidatedRegion> results = Lists.newArrayList();
    private final List<FittedCopyNumber> copyNumbers;

    public SmoothedRegions(final List<FittedCopyNumber> copyNumbers) {
        this.copyNumbers = copyNumbers;
    }

    public List<ConsolidatedRegion> smooth(final List<ConsolidatedRegion> megaRegions) {

        int minIndex = 0;
        int maxIndex = 0;

        for (int i = 0; i < megaRegions.size(); i++) {
            ConsolidatedRegion currentRegion = megaRegions.get(i);
            @Nullable
            ConsolidatedRegion nextRegion = i < megaRegions.size() - 1 ? megaRegions.get(i + 1) : null;

            minIndex = Math.max(minIndex, findChromosomeStartIndex(minIndex, currentRegion.chromosome()));

            int chromosomeStartIndex = findChromosomeStartIndex(minIndex, currentRegion.chromosome());
            minIndex = Math.max(minIndex, chromosomeStartIndex);

            int regionStartIndex = findRegionStartIndex(minIndex, currentRegion);
            int regionEndIndex = findRegionEndIndex(regionStartIndex, currentRegion);

            maxIndex = nextRegion == null ? copyNumbers.size() : findRegionStartIndex(regionEndIndex, nextRegion) - 1;
            minIndex = startConsolidatedRegion(minIndex, regionStartIndex, regionEndIndex, maxIndex) + 1;
        }

        return results;
    }

    private int findChromosomeStartIndex(int minIndex, String chromosome) {
        return indexOf(minIndex, copyNumber -> copyNumber.chromosome().equals(chromosome));
    }

    private int findRegionStartIndex(int minIndex, ConsolidatedRegion region) {
        return indexOf(minIndex, copyNumber -> copyNumber.chromosome().equals(region.chromosome())
                && copyNumber.start() == region.start());
    }

    private int findRegionEndIndex(int minIndex, final ConsolidatedRegion region) {
        return indexOf(minIndex,
                copyNumber -> copyNumber.chromosome().equals(region.chromosome()) && copyNumber.end() == region.end());

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

    private int startConsolidatedRegion(final int minIndex, final int regionStartIndex, final int regionEndIndex,
            final int maxIndex) {
        assert (maxIndex >= regionEndIndex);
        assert (minIndex <= regionStartIndex);
        ConsolidatedRegionBuilder forwardBuilder = new ConsolidatedRegionBuilder(copyNumbers.get(regionStartIndex));

        // Pre-region
        ConsolidatedRegionBuilder reverseBuilder = forwardBuilder;
        final Deque<ConsolidatedRegion> preRegionResults = new ArrayDeque<>();
        for (int i = regionStartIndex; i >= minIndex; i--) {

            FittedCopyNumber copyNumber = copyNumbers.get(i);
            if (isSimilar(copyNumber, reverseBuilder)) {
                reverseBuilder.extendRegion(copyNumber);
            } else {
                if (reverseBuilder != forwardBuilder) {
                    preRegionResults.addFirst(reverseBuilder.build());
                }
                reverseBuilder = new ConsolidatedRegionBuilder(copyNumber);
            }
        }

        // Finalise pre-region
        if (reverseBuilder != forwardBuilder) {
            preRegionResults.addFirst(reverseBuilder.build());
            results.addAll(preRegionResults);
        }

        // Inside region
        for (int i = regionStartIndex + 1; i <= regionEndIndex; i++) {
            FittedCopyNumber copyNumber = copyNumbers.get(i);
            if (isSimilar(copyNumber, forwardBuilder)) {
                forwardBuilder.extendRegion(copyNumber);
            } else {
                results.add(forwardBuilder.build());
                forwardBuilder = new ConsolidatedRegionBuilder(copyNumber);
            }
        }

        // Post-region
        for (int i = regionEndIndex + 1; i <= maxIndex; i++) {
            FittedCopyNumber copyNumber = copyNumbers.get(i);
            if (isSimilar(copyNumber, forwardBuilder)) {
                forwardBuilder.extendRegion(copyNumber);
            } else {
                results.add(forwardBuilder.build());
                return i - 1;
            }
        }

        results.add(forwardBuilder.build());
        return maxIndex;
    }

    private boolean isDiploid(FittedCopyNumber copyNumber) {
        return Doubles.greaterOrEqual(copyNumber.normalCNVRatio(), DIPLOID_MIN_RATIO) && Doubles.lessOrEqual(
                copyNumber.normalCNVRatio(), DIPLOID_MAX_RATIO);
    }


    private boolean isSimilar(FittedCopyNumber copyNumber, ConsolidatedRegionBuilder builder) {

        if (!builder.chromosome().equals(copyNumber.chromosome())) {
            return false;
        }

        if (!isDiploid(copyNumber)) {
            return true;
        }

        if (copyNumber.bafCount() == 0) {
            return true;
        }

        double ratioDevation = Math.abs(copyNumber.ratioOfRatios() - builder.averageRatioOfRatios());
        if (Doubles.greaterThan(ratioDevation, RATIO_RANGE)) {
            return false;
        }

        double bafDeviation = Math.abs(copyNumber.actualBAF() - builder.averageBAF());
        if (Doubles.greaterThan(bafDeviation, allowedBAFDeviation(copyNumber.bafCount()))) {
            return false;
        }

        return true;
    }

    private double allowedBAFDeviation(int bafCount) {
        if (bafCount >= MAX_BAF_COUNT) {
            return MAX_BAF_RANGE;
        }

        return (MIN_BAF_RANGE - MAX_BAF_RANGE) / MAX_BAF_COUNT * bafCount + MAX_BAF_RANGE;
    }

}
