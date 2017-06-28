package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class BaseSmoothedRegions {

    private static final double DIPLOID_MIN_RATIO = 0.75;
    private static final double DIPLOID_MAX_RATIO = 1.25;

    static boolean isDiploid(@NotNull final FittedRegion copyNumber) {
        return Doubles.greaterOrEqual(copyNumber.observedNormalRatio(), DIPLOID_MIN_RATIO) && Doubles.lessOrEqual(
                copyNumber.observedNormalRatio(), DIPLOID_MAX_RATIO);
    }
}
