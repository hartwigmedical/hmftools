package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class LowConfidenceSmoothedRegions {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final List<FittedRegion> fittedRegions;
    @NotNull
    private final List<FittedRegion> smoothedRegions = Lists.newArrayList();
    private final CopyNumberDeviation deviation;

    LowConfidenceSmoothedRegions(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<FittedRegion> fittedRegions) {
        this.purityAdjuster = purityAdjuster;
        this.fittedRegions = fittedRegions;
        this.deviation = new CopyNumberDeviation(purityAdjuster);
        run();
    }

    List<FittedRegion> smoothedRegions() {
        return smoothedRegions;
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
                    smoothedRegions.add(builder.region());
                    builder = new CombinedFittedRegion(false, fittedRegion);
                }
            }
        }

        if (builder != null) {
            smoothedRegions.add(builder.region());
        }
    }

    private boolean isSimilar(@NotNull final FittedRegion newRegion, @NotNull final FittedRegion currentRegion) {
        return newRegion.status().equals(FreecStatus.GERMLINE) || deviation.withinCopyNumberTolerance(currentRegion, newRegion)
                || Doubles.isZero(newRegion.tumorCopyNumber()) || Doubles.isZero(currentRegion.tumorCopyNumber());
    }

}
