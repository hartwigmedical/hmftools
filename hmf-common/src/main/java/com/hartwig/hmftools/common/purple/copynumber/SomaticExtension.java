package com.hartwig.hmftools.common.purple.copynumber;

import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.base.MoreObjects;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class SomaticExtension {

    @NotNull
    static List<FittedRegion> fittedRegions(@NotNull final PurityAdjuster adjuster, @NotNull final Collection<FittedRegion> fittedRegions) {
        final SomaticExtension extension = new SomaticExtension(adjuster, fittedRegions);
        return extension.fittedRegions();
    }

    private static final int MIN_BAF_COUNT_TO_WEIGH_WITH_BAF = 50;

    private static final Logger LOGGER = LogManager.getLogger(SomaticExtension.class);
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    private List<CombinedFittedRegion> regions = Lists.newLinkedList();
    private final BAFDeviation bafDeviation;
    private final CopyNumberDeviation copyNumberDeviation;

    private SomaticExtension(@NotNull final PurityAdjuster adjuster, @NotNull final Collection<FittedRegion> fittedRegions) {
        this.bafDeviation = new BAFDeviation();
        this.copyNumberDeviation = new CopyNumberDeviation(adjuster);
        somaticExtension(fittedRegions);
    }

    @NotNull
    private List<FittedRegion> fittedRegions() {
        return regions.stream().map(CombinedFittedRegion::region).collect(Collectors.toList());
    }

    private void somaticExtension(@NotNull final Collection<FittedRegion> fittedRegions) {
        final boolean bafWeighted = fittedRegions.stream().anyMatch(x -> x.bafCount() >= MIN_BAF_COUNT_TO_WEIGH_WITH_BAF);

        for (FittedRegion fittedRegion : fittedRegions) {
            regions.add(new CombinedFittedRegion(bafWeighted, fittedRegion, false));
        }

        int highestConfidenceIndex = nextIndex();
        while (highestConfidenceIndex > -1) {
            final CombinedFittedRegion highestConfidence = regions.get(highestConfidenceIndex);

            LOGGER.debug("Selected region {}", toString(highestConfidence.region()));
            extendRight(highestConfidenceIndex);
            extendLeft(highestConfidenceIndex);

            LOGGER.debug("Completed region {}", toString(highestConfidence.region()));
            highestConfidence.setModified();
            highestConfidenceIndex = nextIndex();
        }

    }

    private void extendRight(int startIndex) {
        assert (startIndex < regions.size());
        final CombinedFittedRegion target = regions.get(startIndex);
        int targetIndex = startIndex + 1;

        while (targetIndex < regions.size()) {
            final FittedRegion neighbour = regions.get(targetIndex).region();
            if (!merge(target, neighbour)) {
                return;
            }

            regions.remove(targetIndex);
            LOGGER.debug("Merged in right region {}", toString(neighbour));
        }
    }

    private void extendLeft(int startIndex) {
        assert (startIndex < regions.size());
        final CombinedFittedRegion target = regions.get(startIndex);

        int targetIndex = startIndex - 1;
        while (targetIndex >= 0) {
            final FittedRegion neighbour = regions.get(targetIndex).region();
            if (!merge(target, neighbour)) {
                return;
            }

            regions.remove(targetIndex);
            LOGGER.debug("Merged in left region {}", toString(neighbour));
            targetIndex--;
        }

    }

    private boolean breakForCentromereStart(@NotNull final CombinedFittedRegion target, @NotNull final FittedRegion neighbour) {
        if (target.region().start() < neighbour.start()) {
            return neighbour.support() == SegmentSupport.CENTROMERE;
        }

        return target.region().support() == SegmentSupport.CENTROMERE;
    }

    private boolean breakForStructuralVariant(@NotNull final CombinedFittedRegion target, @NotNull final FittedRegion neighbour) {
        if (target.region().start() < neighbour.start()) {
            return neighbour.support() != SegmentSupport.NONE;
        }

        return target.region().support() != SegmentSupport.NONE;
    }

    private boolean merge(@NotNull final CombinedFittedRegion target, @NotNull final FittedRegion neighbour) {

        if (breakForCentromereStart(target, neighbour) || breakForStructuralVariant(target, neighbour)) {
            return false;
        }

        final boolean isNeighbourValid = isValidSomatic(neighbour);

        if (!isNeighbourValid) {
            target.combine(neighbour, false);
            return true;
        } else if (inTolerance(target.region(), neighbour)) {
            target.combine(neighbour, true);
            return true;
        }

        return false;
    }

    private boolean isValidSomatic(@NotNull final FittedRegion region) {
        return region.status() == ObservedRegionStatus.SOMATIC && (region.support() != SegmentSupport.NONE
                || region.observedTumorRatioCount() > 5);
    }

    private boolean inTolerance(@NotNull final FittedRegion left, @NotNull final FittedRegion right) {
        return bafDeviation.inTolerance(left, right) && copyNumberDeviation.inTolerance(left, right);
    }

    private int nextIndex() {

        int indexOfLargestBaf = -1;
        int indexOfLargestLength = -1;

        int largestBAFCount = 0;
        long largestBases = 0;

        for (int i = 0; i < regions.size(); i++) {
            final CombinedFittedRegion combined = regions.get(i);
            final FittedRegion region = combined.region();
            if (!combined.isModified() && region.status().equals(ObservedRegionStatus.SOMATIC)) {

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
                .add("support", region.support())
                .add("copyNumber", FORMAT.format(region.tumorCopyNumber()))
                .toString();
    }
}
