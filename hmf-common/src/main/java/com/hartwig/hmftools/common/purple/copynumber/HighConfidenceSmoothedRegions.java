package com.hartwig.hmftools.common.purple.copynumber;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

class HighConfidenceSmoothedRegions {

    private final List<CombinedFittedRegion> combinedRegions = Lists.newArrayList();
    private final List<? extends GenomeRegion> highConfidenceRegions;
    private final List<FittedRegion> fittedRegions;
    private final CopyNumberDeviation deviation;

    HighConfidenceSmoothedRegions(@NotNull final PurityAdjuster purityAdjuster,
            @NotNull final List<? extends GenomeRegion> highConfidenceRegions, @NotNull final List<FittedRegion> fittedRegions) {
        this.deviation = new CopyNumberDeviation(purityAdjuster);
        this.highConfidenceRegions = highConfidenceRegions;
        this.fittedRegions = fittedRegions;

        run();
    }

    @NotNull
    List<FittedRegion> smoothedRegions() {
        return secondPass(combinedRegions).stream().map(CombinedFittedRegion::region).collect(Collectors.toList());
    }

    private void run() {
        if (!highConfidenceRegions.isEmpty()) {
            int largestIncludedIndex = -1;
            CombinedFittedRegion currentBuilder;

            for (int i = 0; i < highConfidenceRegions.size(); i++) {
                final GenomeRegion currentRegion = highConfidenceRegions.get(i);
                int startOfRegionIndex = indexOfStart(largestIncludedIndex + 1, currentRegion);
                int endOfRegionIndex = indexOfEnd(startOfRegionIndex, currentRegion);

                // JOBA: Start new builder
                currentBuilder = new CombinedFittedRegion(true, fittedRegions.get(startOfRegionIndex));

                // JOBA: Go backwards to previous end
                currentBuilder = backwards(startOfRegionIndex - 1, largestIncludedIndex + 1, currentBuilder);

                // JOBA: Go forwards to end of region
                currentBuilder = forwards(startOfRegionIndex + 1, endOfRegionIndex, currentBuilder);

                boolean isLastBroadRegion = i == highConfidenceRegions.size() - 1;
                if (isLastBroadRegion) {
                    currentBuilder = forwards(endOfRegionIndex + 1, fittedRegions.size() - 1, currentBuilder);
                    combinedRegions.add(currentBuilder);
                } else {
                    int nextStartIndex = indexOfStart(endOfRegionIndex + 1, highConfidenceRegions.get(i + 1));
                    largestIncludedIndex = forwardsUntilDifferent(endOfRegionIndex + 1, nextStartIndex - 1, currentBuilder);
                    combinedRegions.add(currentBuilder);
                }
            }
        }
    }

    private int forwardsUntilDifferent(int startIndex, int endIndex, @NotNull CombinedFittedRegion builder) {

        final List<FittedRegion> buffer = Lists.newArrayList();

        for (int i = startIndex; i <= endIndex; i++) {
            FittedRegion region = fittedRegions.get(i);
            boolean isGermline = isGermline(region);
            boolean isSVSupported = region.structuralVariantSupport() != StructuralVariantSupport.NONE;
            if (isGermline && (isSVSupported || !buffer.isEmpty())) {
                buffer.add(region);
            } else {
                boolean isSimilar = isSimilar(region, builder.region());
                if (isSimilar) {
                    // Flush buffer
                    Iterator<FittedRegion> iterator = buffer.iterator();
                    while (iterator.hasNext()) {
                        FittedRegion next = iterator.next();
                        builder.combine(next);
                        iterator.remove();
                    }

                    builder.combine(region);
                } else {
                    return i - 1 - buffer.size();
                }
            }
        }

        return endIndex - buffer.size();
    }

    @NotNull
    private CombinedFittedRegion forwards(int startIndex, int endIndex, final @NotNull CombinedFittedRegion builder) {
        CombinedFittedRegion current = builder;
        for (int i = startIndex; i <= endIndex; i++) {
            FittedRegion copyNumber = fittedRegions.get(i);
            if (isSimilar(copyNumber, current.region())) {
                current.combine(copyNumber);
            } else {
                combinedRegions.add(current);
                current = new CombinedFittedRegion(true, copyNumber);
            }
        }

        return current;
    }

    @NotNull
    private CombinedFittedRegion backwards(int startIndex, int endIndex, @NotNull final CombinedFittedRegion forwardBuilder) {
        final Deque<CombinedFittedRegion> preRegions = new ArrayDeque<>();
        CombinedFittedRegion reverseBuilder = forwardBuilder;

        for (int i = startIndex; i >= endIndex; i--) {
            final FittedRegion copyNumber = fittedRegions.get(i);
            if (isSimilar(copyNumber, reverseBuilder.region())) {
                reverseBuilder.combine(copyNumber);
            } else {
                if (reverseBuilder != forwardBuilder) {
                    preRegions.addFirst(reverseBuilder);
                }
                reverseBuilder = new CombinedFittedRegion(true, copyNumber);
            }
        }

        if (reverseBuilder != forwardBuilder) {
            preRegions.addFirst(reverseBuilder);
        }

        combinedRegions.addAll(preRegions);
        return forwardBuilder;
    }

    private boolean isGermline(@NotNull final FittedRegion newRegion) {
        return !newRegion.status().equals(FreecStatus.SOMATIC);
    }

    private boolean isSimilar(@NotNull final FittedRegion newRegion, @NotNull final FittedRegion combinedRegion) {
        int bafCount = Math.min(newRegion.bafCount(), combinedRegion.bafCount());
        if (isGermline(newRegion)) {
            return true;
        }

        if (!deviation.withinCopyNumberTolerance(combinedRegion, newRegion)) {
            return false;
        }

        if (bafCount > 0 && !Doubles.isZero(combinedRegion.observedBAF())) {
            double bafDeviation = Math.abs(newRegion.observedBAF() - combinedRegion.observedBAF());
            if (Doubles.greaterThan(bafDeviation, allowedBAFDeviation(bafCount))) {
                return false;
            }
        }

        return true;
    }

    static double allowedBAFDeviation(int bafCount) {
        return 1d / Math.max(1, bafCount) * 0.35 + 0.03;
    }

    private int indexOfEnd(int minIndex, @NotNull GenomeRegion region) {
        return indexOf(minIndex, copyNumber -> copyNumber.end() == region.end());
    }

    private int indexOfStart(int minIndex, @NotNull GenomeRegion region) {
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

    private List<CombinedFittedRegion> secondPass(List<CombinedFittedRegion> regions) {

        final List<CombinedFittedRegion> result = Lists.newArrayList();

        CombinedFittedRegion current = regions.get(0);
        for (int i = 1; i < regions.size(); i++) {
            CombinedFittedRegion next = regions.get(i);
            if (isSimilar(next.region(), current.region())) {
                current.combine(next.region());
            } else {
                result.add(current);
                current = next;
            }
        }
        result.add(current);

        return result;
    }

}
