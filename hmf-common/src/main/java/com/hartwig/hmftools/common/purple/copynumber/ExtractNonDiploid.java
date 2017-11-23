package com.hartwig.hmftools.common.purple.copynumber;

import java.util.Collections;
import java.util.EnumSet;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;

import org.jetbrains.annotations.NotNull;

public class ExtractNonDiploid {

    private static final EnumSet<ObservedRegionStatus> ELIGIBLE = EnumSet.of(ObservedRegionStatus.GERMLINE_AMPLIFICATION,
            ObservedRegionStatus.GERMLINE_HET_DELETION,
            ObservedRegionStatus.GERMLINE_HOM_DELETION);

    @NotNull
    static List<CombinedRegion> extractNonDiploid(@NotNull final List<CombinedRegion> regions) {

        final List<CombinedRegion> result = Lists.newArrayList();
        for (CombinedRegion region : regions) {
            result.addAll(extractNonDiploid(region));
        }

        return result;
    }

    @NotNull
    static List<CombinedRegion> extractNonDiploid(@NotNull final CombinedRegion parent) {
        final List<CombinedRegion> children = Lists.newArrayList();

        double copyNumber = parent.tumorCopyNumber();
        for (final FittedRegion fittedRegion : parent.regions()) {
            if (ELIGIBLE.contains(fittedRegion.status())) {
                final double upperBound = upperBound(fittedRegion);
                final double lowerBound = lowerBound(fittedRegion);

                if (Doubles.lessThan(upperBound, Math.min(0.5, copyNumber))) {
                    children.add(create(fittedRegion, upperBound));
                } else if (Doubles.greaterThan(lowerBound, copyNumber)) {
                    children.add(create(fittedRegion, lowerBound));
                }
            }
        }

        return selleysNoMoreGaps(parent, children);
    }

    static CombinedRegion create(@NotNull final FittedRegion child, double newCopyNumber) {
        final CombinedRegion result = new CombinedRegion(false, child, false);
        result.setTumorCopyNumber(CopyNumberMethod.NON_DIPLOID, newCopyNumber);
        return result;
    }

    @NotNull
    static List<CombinedRegion> selleysNoMoreGaps(@NotNull final CombinedRegion parent, final List<CombinedRegion> children) {
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

    private static double upperBound(@NotNull final FittedRegion region) {
        return Math.max(region.tumorCopyNumber(), region.refNormalisedCopyNumber()) / (2 * Math.min(1, region.observedNormalRatio()));
    }

    private static double lowerBound(@NotNull final FittedRegion region) {
        return Math.min(region.tumorCopyNumber(), region.refNormalisedCopyNumber()) / (2 * Math.max(1, region.observedNormalRatio()));
    }
}
