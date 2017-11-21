package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;

import org.jetbrains.annotations.NotNull;

@Deprecated
class LowConfidenceSmoothedRegions {

    @NotNull
    private final List<FittedRegion> fittedRegions;
    @NotNull
    private final List<CombinedRegion> smoothedRegions = Lists.newArrayList();
    private final CopyNumberDeviation deviation;

    LowConfidenceSmoothedRegions(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<FittedRegion> fittedRegions) {
        this.fittedRegions = fittedRegions;
        this.deviation = new CopyNumberDeviation(purityAdjuster);
        run();
    }

    List<FittedRegion> smoothedRegions() {
        return finalPass(smoothedRegions).stream()
                .map(CombinedRegion::region)
                .collect(Collectors.toList());
    }

    private void run() {

        CombinedRegion builder = null;
        for (FittedRegion fittedRegion : fittedRegions) {

            if (builder == null) {
                builder = new CombinedRegion(false, fittedRegion);
            } else {
                if (isSimilar(fittedRegion, builder.region())) {
                    builder.combine(fittedRegion);
                } else {
                    smoothedRegions.add(builder);
                    builder = new CombinedRegion(false, fittedRegion);
                }
            }
        }

        if (builder != null) {
            smoothedRegions.add(builder);
        }
    }

    private boolean isSimilar(@NotNull final FittedRegion newRegion, @NotNull final FittedRegion currentRegion) {
        return !newRegion.status().equals(ObservedRegionStatus.SOMATIC) || deviation.inTolerance(currentRegion, newRegion)
                || Doubles.isZero(newRegion.tumorCopyNumber()) || Doubles.isZero(currentRegion.tumorCopyNumber());
    }

    private List<CombinedRegion> finalPass(List<CombinedRegion> regions) {
        return CombinedFittedRegions.mergeLeft(regions, (left, right) -> isSimilar(right.region(), left.region()));
    }
}
