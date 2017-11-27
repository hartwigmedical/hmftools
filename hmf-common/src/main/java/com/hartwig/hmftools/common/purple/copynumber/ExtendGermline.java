package com.hartwig.hmftools.common.purple.copynumber;

import java.util.Collections;
import java.util.EnumSet;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;

import org.jetbrains.annotations.NotNull;

class ExtendGermline {

    private static final double AMPLIFICATION_TOLERANCE = 1;

    private final Gender gender;

    ExtendGermline(final Gender gender) {
        this.gender = gender;
    }

    @NotNull
    List<CombinedRegion> extendGermlineAmplifications(@NotNull final List<CombinedRegion> regions) {
        final EnumSet<ObservedRegionStatus> eligibleStatus = EnumSet.of(ObservedRegionStatus.GERMLINE_AMPLIFICATION);
        return extendGermline(eligibleStatus, regions);
    }

    @NotNull
    List<CombinedRegion> extendGermlineAmplificationsAndDeletions(@NotNull final List<CombinedRegion> regions) {
        final EnumSet<ObservedRegionStatus> eligibleStatus = EnumSet.of(ObservedRegionStatus.GERMLINE_AMPLIFICATION,
                ObservedRegionStatus.GERMLINE_HET_DELETION,
                ObservedRegionStatus.GERMLINE_HOM_DELETION);
        return extendGermline(eligibleStatus, regions);
    }

    @NotNull
    List<CombinedRegion> extractGermlineDeletions(@NotNull final List<CombinedRegion> regions) {
        final EnumSet<ObservedRegionStatus> eligibleStatus =
                EnumSet.of(ObservedRegionStatus.GERMLINE_HET_DELETION, ObservedRegionStatus.GERMLINE_HOM_DELETION);

        final List<CombinedRegion> result = Lists.newArrayList();
        for (CombinedRegion parent : regions) {
            result.addAll(extractChildren(eligibleStatus, parent));
        }

        return result;
    }

    @NotNull
    private List<CombinedRegion> extendGermline(@NotNull final EnumSet<ObservedRegionStatus> eligbleStatus,
            @NotNull final List<CombinedRegion> regions) {

        final List<CombinedRegion> result = Lists.newArrayList();
        for (CombinedRegion parent : regions) {
            final List<CombinedRegion> children = extractChildren(eligbleStatus, parent);
            if (children.isEmpty()) {
                result.add(parent);
            } else {
                result.addAll(selleysNoMoreGaps(parent, children));
            }
        }
        return result;
    }

    @NotNull
    private List<CombinedRegion> extractChildren(@NotNull final EnumSet<ObservedRegionStatus> eligbleStatus,
            @NotNull final CombinedRegion parent) {
        final List<CombinedRegion> children = Lists.newArrayList();

        double baf = parent.tumorBAF();
        double copyNumber = parent.tumorCopyNumber();
        for (final FittedRegion fittedRegion : parent.regions()) {
            if (eligbleStatus.contains(fittedRegion.status())) {
                if (fittedRegion.status().equals(ObservedRegionStatus.GERMLINE_AMPLIFICATION)) {
                    final double lowerBound = lowerBound(fittedRegion);
                    if (Doubles.greaterThan(lowerBound, copyNumber + AMPLIFICATION_TOLERANCE)) {
                        children.add(createChild(fittedRegion, lowerBound, baf));
                    }
                }

                if (fittedRegion.status().equals(ObservedRegionStatus.GERMLINE_HET_DELETION)) {
                    final double upperBound = upperBound(fittedRegion);
                    if (Doubles.lessThan(upperBound, Math.min(0.5, copyNumber))) {
                        children.add(createChild(fittedRegion, upperBound, baf));
                    }
                }

                if (fittedRegion.status().equals(ObservedRegionStatus.GERMLINE_HOM_DELETION)) {
                    children.add(createChild(fittedRegion, fittedRegion.refNormalisedCopyNumber(), baf));
                }
            }
        }

        return extendRight(children);
    }

    @NotNull
    static CombinedRegion createChild(@NotNull final FittedRegion child, double newCopyNumber, double newBaf) {
        final CombinedRegion result = new CombinedRegion(false, child, false);
        result.setTumorCopyNumber(method(child), newCopyNumber);
        result.setInferredTumorBAF(newBaf);
        return result;
    }

    @NotNull
    private static CopyNumberMethod method(@NotNull final FittedRegion child) {
        switch (child.status()) {
            case GERMLINE_HET_DELETION:
                return CopyNumberMethod.GERMLINE_HET2HOM_DELETION;
            case GERMLINE_HOM_DELETION:
                return CopyNumberMethod.GERMLINE_HOM_DELETION;
            default:
                return CopyNumberMethod.GERMLINE_AMPLIFICATION;
        }
    }

    @NotNull
    private static List<CombinedRegion> selleysNoMoreGaps(@NotNull final CombinedRegion parent, final List<CombinedRegion> children) {
        if (children.isEmpty()) {
            return Collections.singletonList(parent);
        }

        final List<CombinedRegion> result = Lists.newArrayList();
        long nextStart = parent.start();
        for (CombinedRegion child : children) {
            if (child.start() > nextStart) {
                result.add(CombinedRegionFactory.reduce(parent, nextStart, child.start() - 1));
            }

            result.add(child);
            nextStart = child.end() + 1;
        }

        if (nextStart <= parent.end()) {
            result.add(CombinedRegionFactory.reduce(parent, nextStart, parent.end()));
        }

        return result;
    }

    @NotNull
    private static List<CombinedRegion> extendRight(@NotNull final List<CombinedRegion> children) {

        int i = 0;
        while (i < children.size() - 1) {
            final CombinedRegion target = children.get(i);
            final CombinedRegion neighbour = children.get(i + 1);

            if (target.region().status().equals(neighbour.region().status()) && target.end() + 1 == neighbour.start()) {
                target.extendWithUnweightedAverage(neighbour.region());
                children.remove(i + 1);
            } else {
                i++;
            }
        }

        return children;
    }

    private double upperBound(@NotNull final FittedRegion region) {
        double expectedNormalRatio = HumanChromosome.fromString(region.chromosome()).isHomologous(gender) ? 1 : 0.5;
        return Math.max(region.tumorCopyNumber(), region.refNormalisedCopyNumber()) / (2 * Math.min(expectedNormalRatio,
                region.observedNormalRatio()));
    }

    private double lowerBound(@NotNull final FittedRegion region) {
        double expectedNormalRatio = HumanChromosome.fromString(region.chromosome()).isHomologous(gender) ? 1 : 0.5;
        return Math.min(region.tumorCopyNumber(), region.refNormalisedCopyNumber()) / (2 * Math.max(expectedNormalRatio,
                region.observedNormalRatio()));
    }
}
