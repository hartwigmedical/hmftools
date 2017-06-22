package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class LowConfidenceSmoothedRegions extends BaseSmoothedRegions {

    private final double purity;
    private final List<FittedRegion> fittedRegions;
    private final List<PurpleCopyNumber> smoothedRegions = Lists.newArrayList();

    LowConfidenceSmoothedRegions(final double purity, final List<FittedRegion> fittedRegions) {
        this.purity = purity;
        this.fittedRegions = fittedRegions;
        run();
    }

    List<PurpleCopyNumber> smoothedRegions() {
        return smoothedRegions;
    }

    private void run() {

        LowConfidenceCopyNumberBuilder builder = null;
        for (FittedRegion fittedRegion : fittedRegions) {

            if (builder == null) {
                builder = new LowConfidenceCopyNumberBuilder(purity, fittedRegion);
            } else {
                if (isSimilar(fittedRegion, builder)) {
                    builder.extendRegion(fittedRegion);
                } else {
                    smoothedRegions.add(builder.build());
                    builder = new LowConfidenceCopyNumberBuilder(purity, fittedRegion);
                }
            }
        }

        if (builder != null) {
            smoothedRegions.add(builder.build());
        }
    }

    private static boolean isSimilar(@NotNull final FittedRegion copyNumber, @NotNull final LowConfidenceCopyNumberBuilder builder) {
        int bafCount = copyNumber.bafCount();
        if (!isDiploid(copyNumber)) {
            return true;
        }

        double tumorCopyNumberDeviation = Math.abs(copyNumber.tumorCopyNumber() - builder.averageTumorCopyNumber());
        double refNormalisedCopyNumberDeviation = Math.abs(copyNumber.refNormalisedCopyNumber() - builder.averageRefNormalisedCopyNumber());
        double copyNumberDeviation = Math.min(tumorCopyNumberDeviation, refNormalisedCopyNumberDeviation);
        if (Doubles.greaterThan(copyNumberDeviation, allowedCopyNumberDeviation(bafCount))) {
            return false;
        }

        return true;
    }

}
