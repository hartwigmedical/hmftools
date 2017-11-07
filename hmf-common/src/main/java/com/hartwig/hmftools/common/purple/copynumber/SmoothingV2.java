package com.hartwig.hmftools.common.purple.copynumber;

import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.base.MoreObjects;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class SmoothingV2 {
    private static final int MIN_BAF_COUNT_TO_WEIGH_WITH_BAF = 50;

    private static final Logger LOGGER = LogManager.getLogger(SmoothingV2.class);
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    private List<CombinedFittedRegion> regions = Lists.newLinkedList();
    private final BAFDeviation bafDeviation;
    private final CopyNumberDeviation copyNumberDeviation;

    SmoothingV2(@NotNull final PurityAdjuster adjuster) {
        this.bafDeviation = new BAFDeviation();
        this.copyNumberDeviation = new CopyNumberDeviation(adjuster);
    }

    @NotNull
    private List<FittedRegion> smoothedRegions() {
        return regions.stream().map(CombinedFittedRegion::region).collect(Collectors.toList());
    }

    @NotNull
    List<FittedRegion> smooth(@NotNull final Collection<FittedRegion> fittedRegions) {

        final boolean bafWeighted = fittedRegions.stream().anyMatch(x -> x.bafCount() >= MIN_BAF_COUNT_TO_WEIGH_WITH_BAF);

        for (FittedRegion fittedRegion : fittedRegions) {
            regions.add(new CombinedFittedRegion(bafWeighted, fittedRegion, false));
        }

        int highestConfidenceIndex = nextIndex();
        while (highestConfidenceIndex > -1) {
            final CombinedFittedRegion highestConfidence = regions.get(highestConfidenceIndex);

            LOGGER.info("Selected region {}", toString(highestConfidence.region()));
            smoothRight(highestConfidenceIndex);
            smoothLeft(highestConfidenceIndex);

            LOGGER.info("Completed region {}", toString(highestConfidence.region()));
            highestConfidence.setModified();
            highestConfidenceIndex = nextIndex();
        }

        return smoothedRegions();
    }

    private void smoothRight(int startIndex) {
        assert (startIndex < regions.size());
        final CombinedFittedRegion target = regions.get(startIndex);
        int targetIndex = startIndex + 1;

        while (targetIndex < regions.size()) {
            final FittedRegion neighbour = regions.get(targetIndex).region();
            initialiseTarget(target, neighbour);
            if (!mergeIfPossible(target, neighbour)) {
                return;
            }

            regions.remove(targetIndex);
            LOGGER.info("Merged in right region {}", toString(neighbour));
        }
    }

    private void smoothLeft(int startIndex) {
        assert (startIndex < regions.size());
        final CombinedFittedRegion target = regions.get(startIndex);

        int targetIndex = startIndex - 1;
        while (targetIndex >= 0) {
            final FittedRegion neighbour = regions.get(targetIndex).region();
            initialiseTarget(target, neighbour);

            if (!mergeIfPossible(target, neighbour)) {
                return;
            }

            regions.remove(targetIndex);
            LOGGER.info("Merged in left region {}", toString(neighbour));
            targetIndex--;
        }

    }

    private void initialiseTarget(@NotNull final CombinedFittedRegion target, @NotNull final FittedRegion neighbour) {
        final FittedRegion targetRegion = target.region();
        final double neighbourCopyNumber = neighbour.tumorCopyNumber();
        if (!target.isModified() && !(isValidSomatic(targetRegion) || isValidGermline(neighbourCopyNumber, targetRegion))) {
            target.clearValues();
        }
    }

    private boolean breakForCentromereStart(@NotNull final CombinedFittedRegion target, @NotNull final FittedRegion neighbour) {
        if (target.region().start() < neighbour.start()) {
            return neighbour.status() == ObservedRegionStatus.CENTROMERE;
        }

        return target.spansCentromere();
    }

    private boolean mergeIfPossible(@NotNull final CombinedFittedRegion target, @NotNull final FittedRegion neighbour) {

        if (breakForCentromereStart(target, neighbour)) {
            return false;
        }

        final double targetCopyNumber = target.region().tumorCopyNumber();
        final boolean isNeighbourValid = isValidSomatic(neighbour) || isValidGermline(targetCopyNumber, neighbour);
        if (!isNeighbourValid) {
            target.combine(neighbour, false);
        } else if (inTolerance(target.region(), neighbour)) {
            target.combine(neighbour, true);
        } else {
            return false;
        }

        return true;
    }

    private boolean isValidSomatic(@NotNull final FittedRegion region) {
        return region.status() == ObservedRegionStatus.SOMATIC && (region.structuralVariantSupport() != StructuralVariantSupport.NONE
                || region.observedTumorRatioCount() > 5);
    }

    private boolean isValidGermline(double neighbourCopyNumber, @NotNull final FittedRegion region) {
        return isValidGermlineHetrozygous(region) || isValidGermlineAmplification(neighbourCopyNumber, region);
    }

    private boolean isValidGermlineHetrozygous(@NotNull final FittedRegion region) {
        return false;
        //        return region.status() == ObservedRegionStatus.GERMLINE_HET_DELETION && Doubles.lessThan(region.refNormalisedCopyNumber(), 0.5);
    }

    private boolean isValidGermlineAmplification(double neighbourCopyNumber, @NotNull final FittedRegion region) {
        return false;
        //        double tumorCopyNumber = region.refNormalisedCopyNumber() / (Math.ceil(2 * region.observedNormalRatio()));
        //        return region.status() == ObservedRegionStatus.GERMLINE_AMPLIFICATION && Doubles.greaterThan(tumorCopyNumber,
        //                surroundingTumorCopyNumber);
    }

    private boolean inTolerance(@NotNull final FittedRegion target, @NotNull final FittedRegion neighbour) {
        final FittedRegion right = neighbour.start() > target.start() ? neighbour : target;
        return right.structuralVariantSupport() == StructuralVariantSupport.NONE
                && inTolerances(target, neighbour);
    }

    private boolean inTolerances(@NotNull final FittedRegion left, @NotNull final FittedRegion right) {
        return bafDeviation.inTolerance(left, right) && copyNumberDeviation.inTolerance(left, right);
    }

    private int nextIndex() {
        int nextSomatic = nextIndex(x -> x == ObservedRegionStatus.SOMATIC);
        return nextSomatic == -1 ? nextIndex(x -> true) : nextSomatic;
    }

    private int nextIndex(@NotNull final Predicate<ObservedRegionStatus> statusPredicate) {

        int indexOfLargestBaf = -1;
        int indexOfLargestLength = -1;

        int largestBAFCount = 0;
        long largestBases = 0;

        for (int i = 0; i < regions.size(); i++) {
            final CombinedFittedRegion combined = regions.get(i);
            final FittedRegion region = combined.region();
            if (!combined.isModified() && statusPredicate.test(region.status())) {

                if (region.bafCount() > largestBAFCount) {
                    largestBAFCount = region.bafCount();
                    indexOfLargestBaf = i;
                }

                if (region.bases() > largestBases) {
                    largestBases = region.bases();
                    indexOfLargestLength = i;
                }
            }
        }

        return indexOfLargestBaf > -1 ? indexOfLargestBaf : indexOfLargestLength;
    }

    private String toString(FittedRegion region) {
        return MoreObjects.toStringHelper("FittedRegion")
                .omitNullValues()
                .add("chromosome", region.chromosome())
                .add("start", region.start())
                .add("end", region.end())
                .add("status", region.status())
                .add("sv", region.structuralVariantSupport())
                .add("copyNumber", FORMAT.format(region.tumorCopyNumber()))
                .toString();
    }
}
