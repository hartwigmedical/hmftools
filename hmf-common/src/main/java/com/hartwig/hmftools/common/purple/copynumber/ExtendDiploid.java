package com.hartwig.hmftools.common.purple.copynumber;

import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.function.IntUnaryOperator;

import com.google.common.base.MoreObjects;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class ExtendDiploid {

    private enum Direction {
        LEFT(index -> index - 1),
        RIGHT(index -> index + 1);

        private final IntUnaryOperator indexOperator;

        Direction(final IntUnaryOperator indexOperator) {
            this.indexOperator = indexOperator;
        }

        int moveIndex(int index) {
            return indexOperator.applyAsInt(index);
        }

    }

    private static final int MIN_BAF_COUNT_TO_WEIGH_WITH_BAF = 50;
    private static final Logger LOGGER = LogManager.getLogger(ExtendDiploid.class);
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    private final int minTumorCount;
    private final int minTumorCountAtCentromere;
    private final BAFDeviation bafDeviation;
    private final CopyNumberDeviation copyNumberDeviation;

    ExtendDiploid(@NotNull final PurityAdjuster adjuster, final int minTumorCount, final int minTumorCountAtCentromere) {
        this.minTumorCount = minTumorCount;
        this.minTumorCountAtCentromere = minTumorCountAtCentromere;
        this.bafDeviation = new BAFDeviation();
        this.copyNumberDeviation = new CopyNumberDeviation(adjuster);
    }

    @NotNull
    List<CombinedRegion> extendDiploid(@NotNull final Collection<FittedRegion> fittedRegions) {
        final boolean bafWeighted = fittedRegions.stream().anyMatch(x -> x.bafCount() >= MIN_BAF_COUNT_TO_WEIGH_WITH_BAF);

        final List<CombinedRegion> regions = Lists.newLinkedList();

        for (FittedRegion fittedRegion : fittedRegions) {
            regions.add(new CombinedRegion(bafWeighted, fittedRegion));
        }

        int highestConfidenceIndex = nextIndex(regions);
        while (highestConfidenceIndex > -1) {
            final CombinedRegion highestConfidence = regions.get(highestConfidenceIndex);
            highestConfidence.setCopyNumberMethod(CopyNumberMethod.BAF_WEIGHTED);

            LOGGER.debug("Selected region {}", toString(highestConfidence.region()));
            extendRight(regions, highestConfidenceIndex);
            extendLeft(regions, highestConfidenceIndex);

            LOGGER.debug("Completed region {}", toString(highestConfidence.region()));
            highestConfidenceIndex = nextIndex(regions);
        }

        return regions;
    }

    private void extendRight(@NotNull final List<CombinedRegion> regions, int targetIndex) {
        assert (targetIndex < regions.size());
        int neighbourIndex = targetIndex + 1;

        while (neighbourIndex < regions.size()) {
            if (!merge(regions, Direction.RIGHT, targetIndex)) {
                return;
            }
            regions.remove(neighbourIndex);
        }
    }

    private void extendLeft(@NotNull final List<CombinedRegion> regions, final int targetIndex) {
        assert (targetIndex < regions.size());
        int neighbourIndex = targetIndex - 1;

        while (neighbourIndex >= 0) {
            if (!merge(regions, Direction.LEFT, neighbourIndex + 1)) {
                return;
            }
            regions.remove(neighbourIndex);
            neighbourIndex--;
        }
    }

    private boolean merge(@NotNull final List<CombinedRegion> regions, @NotNull final Direction direction, int targetIndex) {

        final CombinedRegion target = regions.get(targetIndex);
        final FittedRegion neighbour = regions.get(direction.moveIndex(targetIndex)).region();

        if (Extend.doNotExtend(target, neighbour)) {
            return false;
        }

        final boolean isNeighbourDubious = isDubious(neighbour);
        if (isNeighbourDubious) {
            int minTumorCount = nextBigBreakIsCentromere(regions, direction, targetIndex) ? minTumorCountAtCentromere : this.minTumorCount;
            if (inTolerance(target.region(), neighbour)) {
                target.extendWithBAFWeightedAverage(neighbour);
                return true;
            } else if (pushThroughDubiousRegion(minTumorCount, regions, direction, targetIndex)) {
                target.extend(neighbour);
                return true;
            } else {
                return false;
            }
        }

        final boolean isNeighbourValid = isValid(neighbour);
        if (!isNeighbourValid) {
            target.extend(neighbour);
            return true;
        } else if (inTolerance(target.region(), neighbour)) {
            target.extendWithBAFWeightedAverage(neighbour);
            return true;
        }

        return false;
    }

    private boolean isValid(@NotNull final FittedRegion region) {
        return region.status() == GermlineStatus.DIPLOID && (region.support() != SegmentSupport.NONE
                || region.observedTumorRatioCount() >= minTumorCount);
    }

    private boolean isDubious(@NotNull final FittedRegion region) {
        return region.status() == GermlineStatus.DIPLOID && region.support() == SegmentSupport.NONE
                && region.observedTumorRatioCount() < minTumorCount;
    }

    private boolean nextBigBreakIsCentromere(@NotNull final List<CombinedRegion> regions, @NotNull final Direction direction,
            int targetIndex) {
        for (int i = direction.moveIndex(targetIndex); i >= 0 && i < regions.size(); i = direction.moveIndex(i)) {
            final FittedRegion neighbour = regions.get(i).region();
            if (neighbour.support() == SegmentSupport.CENTROMERE) {
                return true;
            }
            if (neighbour.support().isSV()) {
                return false;
            }
        }

        return false;
    }

    private boolean pushThroughDubiousRegion(int minTumorCount, @NotNull final List<CombinedRegion> regions,
            @NotNull final Direction direction, int targetIndex) {

        int dubiousCount = 0;
        final CombinedRegion target = regions.get(targetIndex);
        for (int i = direction.moveIndex(targetIndex); i >= 0 && i < regions.size(); i = direction.moveIndex(i)) {
            final FittedRegion neighbour = regions.get(i).region();

            if (Extend.breakForCentromereStart(target, neighbour)) {
                return dubiousCount < minTumorCount;
            }

            if (Extend.breakForStructuralVariant(target, neighbour)) {
                return false;
            }

            if (isDubious(neighbour)) {
                dubiousCount += neighbour.observedTumorRatioCount();
            }

            if (dubiousCount >= minTumorCount) {
                return false;
            }

            if (isValid(neighbour)) {
                return inTolerance(target.region(), neighbour);
            }
        }

        return false;
    }

    private boolean inTolerance(@NotNull final FittedRegion left, @NotNull final FittedRegion right) {
        return bafDeviation.inTolerance(left, right) && copyNumberDeviation.inTolerance(left, right);
    }

    private static int nextIndex(@NotNull final List<CombinedRegion> regions) {

        int indexOfLargestBaf = -1;
        int indexOfLargestLength = -1;

        int largestBAFCount = 0;
        long largestBases = 0;

        for (int i = 0; i < regions.size(); i++) {
            final CombinedRegion combined = regions.get(i);
            final FittedRegion region = combined.region();
            if (!combined.isProcessed() && region.status().equals(GermlineStatus.DIPLOID)) {

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

    private static String toString(FittedRegion region) {
        return MoreObjects.toStringHelper("FittedRegion")
                .omitNullValues()
                .add("chromosome", region.chromosome())
                .add("start", region.start())
                .add("end", region.end())
                .add("status", region.status())
                .add("support", region.support())
                .add("copyNumber", FORMAT.format(region.tumorCopyNumber()))
                .toString();
    }
}
