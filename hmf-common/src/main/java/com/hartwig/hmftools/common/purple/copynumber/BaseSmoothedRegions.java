package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class BaseSmoothedRegions {

    private static final double DIPLOID_MIN_RATIO = 0.75;
    private static final double DIPLOID_MAX_RATIO = 1.25;

    private static final double MIN_COPY_NUMBER_RANGE = 0.3;
    private static final double MAX_COPY_NUMBER_RANGE = 1.3;

    static boolean isDiploid(@NotNull final FittedRegion copyNumber) {
        return Doubles.greaterOrEqual(copyNumber.observedNormalRatio(), DIPLOID_MIN_RATIO) && Doubles.lessOrEqual(
                copyNumber.observedNormalRatio(), DIPLOID_MAX_RATIO);
    }

    static double allowedCopyNumberDeviation(int bafCount) {
        if (bafCount >= 10) {
            return MIN_COPY_NUMBER_RANGE;
        }
        return (MIN_COPY_NUMBER_RANGE - MAX_COPY_NUMBER_RANGE) / 10 * bafCount + MAX_COPY_NUMBER_RANGE;
    }


}
