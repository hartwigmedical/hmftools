package com.hartwig.hmftools.common.purple.copynumber;

import java.util.Collections;
import java.util.EnumSet;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ModifiableFittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ExtendGermline {

    private static final double AMPLIFICATION_TOLERANCE = 1;
    private static final double MIN_AMPLIFICATION_COPYNUMBER = 5;

    private final Gender gender;

    ExtendGermline(final Gender gender) {
        this.gender = gender;
    }

    @NotNull
    List<CombinedRegion> extendGermlineAmplifications(@NotNull final List<CombinedRegion> regions) {
        final EnumSet<GermlineStatus> eligibleStatus = EnumSet.of(GermlineStatus.AMPLIFICATION);
        return extendGermline(eligibleStatus, regions);
    }

    @NotNull
    List<CombinedRegion> extractGermlineDeletions(@NotNull final List<CombinedRegion> regions) {
        final EnumSet<GermlineStatus> eligibleStatus = EnumSet.of(GermlineStatus.HET_DELETION, GermlineStatus.HOM_DELETION);

        final List<CombinedRegion> result = Lists.newArrayList();
        for (int i = 0; i < regions.size(); i++) {
            final CombinedRegion parent = regions.get(i);
            final CombinedRegion neighbour = i + 1 == regions.size() ? null : regions.get(i + 1);
            result.addAll(extractChildren(eligibleStatus, parent, neighbour));
        }

        return result;
    }

    @NotNull
    private List<CombinedRegion> extendGermline(@NotNull final EnumSet<GermlineStatus> eligibleStatus,
            @NotNull final List<CombinedRegion> regions) {
        final List<CombinedRegion> result = Lists.newArrayList();
        for (int i = 0; i < regions.size(); i++) {
            final CombinedRegion parent = regions.get(i);
            final CombinedRegion neighbour = i + 1 == regions.size() ? null : regions.get(i + 1);
            final List<CombinedRegion> children = extractChildren(eligibleStatus, parent, neighbour);
            if (children.isEmpty()) {
                result.add(parent);
            } else {
                result.addAll(selleysNoMoreGaps(parent, children));
            }
        }
        return result;
    }

    @NotNull
    private List<CombinedRegion> extractChildren(@NotNull final EnumSet<GermlineStatus> eligibleStatus,
            @NotNull final CombinedRegion parent, @Nullable final CombinedRegion neighbour) {
        final List<CombinedRegion> children = Lists.newArrayList();

        double baf = parent.tumorBAF();
        double copyNumber = parent.tumorCopyNumber();
        for (int i = 0; i < parent.regions().size(); i++) {
            final FittedRegion child = parent.regions().get(i);
            final FittedRegion next = i + 1 == parent.regions().size()
                    ? (neighbour == null ? null : neighbour.regions().get(0))
                    : parent.regions().get(i + 1);

            if (eligibleStatus.contains(child.status())) {
                if (child.status().equals(GermlineStatus.AMPLIFICATION)) {
                    final double lowerBound = lowerBound(child);
                    if (isValidAmplification(copyNumber, lowerBound, child, next)) {
                        children.add(createChild(child, lowerBound, baf));
                    }
                }

                if (child.status().equals(GermlineStatus.HET_DELETION)) {
                    final double upperBound = upperBound(child);
                    if (Doubles.lessThan(upperBound, Math.min(0.5, copyNumber))) {
                        children.add(createChild(child, upperBound, baf));
                    }
                }

                if (child.status().equals(GermlineStatus.HOM_DELETION)) {
                    children.add(createChild(child, child.refNormalisedCopyNumber(), baf));
                }
            }
        }

        return extendRight(children);
    }

    static boolean isValidAmplification(double parentCopyNumber, double lowerBound, @NotNull final FittedRegion child,
            @Nullable final FittedRegion next) {
        boolean adjacentToCentromere =
                child.support() == SegmentSupport.CENTROMERE || (next != null && next.support() == SegmentSupport.CENTROMERE);
        boolean adjacentToSV = child.support().isSV() || (next != null && next.support().isSV());
        return adjacentToSV && !adjacentToCentromere && Doubles.greaterOrEqual(lowerBound, MIN_AMPLIFICATION_COPYNUMBER)
                && Doubles.greaterOrEqual(lowerBound, parentCopyNumber + AMPLIFICATION_TOLERANCE);
    }

    @NotNull
    private static CombinedRegion createChild(@NotNull final FittedRegion child, double newCopyNumber, double newBaf) {
        final CombinedRegion result = new BafWeightedRegion(child);
        result.setTumorCopyNumber(method(child), newCopyNumber);
        result.setInferredTumorBAF(newBaf);
        return result;
    }

    @NotNull
    private static CopyNumberMethod method(@NotNull final FittedRegion child) {
        switch (child.status()) {
            case HET_DELETION:
                return CopyNumberMethod.GERMLINE_HET2HOM_DELETION;
            case HOM_DELETION:
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
                result.add(reduce(parent, nextStart, child.start() - 1));
            }

            result.add(child);
            nextStart = child.end() + 1;
        }

        if (nextStart <= parent.end()) {
            result.add(reduce(parent, nextStart, parent.end()));
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
        double expectedNormalRatio = HumanChromosome.fromString(region.chromosome()).isDiploid(gender) ? 1 : 0.5;
        return Math.max(region.tumorCopyNumber(), region.refNormalisedCopyNumber()) / (2 * Math.min(expectedNormalRatio,
                region.observedNormalRatio()));
    }

    private double lowerBound(@NotNull final FittedRegion region) {
        double expectedNormalRatio = HumanChromosome.fromString(region.chromosome()).isDiploid(gender) ? 1 : 0.5;
        return Math.min(region.tumorCopyNumber(), region.refNormalisedCopyNumber()) / (2 * Math.max(expectedNormalRatio,
                region.observedNormalRatio()));
    }

    @NotNull
    private static CombinedRegion reduce(@NotNull final CombinedRegion parent, long start, long end) {
        assert (start >= parent.start());
        assert (end <= parent.end());

        int bafCount = 0;
        int depthWindowCount = 0;
        for (FittedRegion fittedRegion : parent.regions()) {
            if (fittedRegion.start() >= start && fittedRegion.end() <= end) {
                bafCount += fittedRegion.bafCount();
                depthWindowCount += fittedRegion.depthWindowCount();
            }
        }

        final ModifiableFittedRegion smallerRegion = ModifiableFittedRegion.create()
                .from(parent.region())
                .setStart(start)
                .setEnd(end)
                .setBafCount(bafCount)
                .setDepthWindowCount(depthWindowCount);

        CombinedRegion result = new BafWeightedRegion(smallerRegion);
        result.setCopyNumberMethod(parent.copyNumberMethod());

        for (FittedRegion fittedRegion : parent.regions()) {
            if (fittedRegion.start() >= start && fittedRegion.end() <= end) {
                result.extend(fittedRegion);
            }
        }

        return result;
    }
}
