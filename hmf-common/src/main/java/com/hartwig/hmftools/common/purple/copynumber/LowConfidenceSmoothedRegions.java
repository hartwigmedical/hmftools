package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;

import org.jetbrains.annotations.NotNull;

class LowConfidenceSmoothedRegions {

    @NotNull
    private final List<FittedRegion> fittedRegions;
    @NotNull
    private final List<CombinedFittedRegion> smoothedRegions = Lists.newArrayList();
    private final CopyNumberDeviation deviation;

    LowConfidenceSmoothedRegions(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<FittedRegion> fittedRegions) {
        this.fittedRegions = fittedRegions;
        this.deviation = new CopyNumberDeviation(purityAdjuster);
        run();
    }

    List<FittedRegion> smoothedRegions() {
        return finalPass(smoothedRegions).stream()
                .map(CombinedFittedRegion::region)
                .collect(Collectors.toList());
    }

    private void run() {

        CombinedFittedRegion builder = null;
        for (FittedRegion fittedRegion : fittedRegions) {

            if (builder == null) {
                builder = new CombinedFittedRegion(false, fittedRegion);
            } else {
                if (isSimilar(fittedRegion, builder.region())) {
                    builder.combine(fittedRegion);
                } else {
                    smoothedRegions.add(builder);
                    builder = new CombinedFittedRegion(false, fittedRegion);
                }
            }
        }

        if (builder != null) {
            smoothedRegions.add(builder);
        }
    }

    private boolean isSimilar(@NotNull final FittedRegion newRegion, @NotNull final FittedRegion currentRegion) {
        return newRegion.status().equals(ObservedRegionStatus.GERMLINE) || deviation.withinCopyNumberTolerance(currentRegion, newRegion)
                || Doubles.isZero(newRegion.tumorCopyNumber()) || Doubles.isZero(currentRegion.tumorCopyNumber());
    }

    private List<CombinedFittedRegion> finalPass(List<CombinedFittedRegion> regions) {
        return CombinedFittedRegions.mergeLeft(regions, (left, right) -> isSimilar(right.region(), left.region()));
    }
}
