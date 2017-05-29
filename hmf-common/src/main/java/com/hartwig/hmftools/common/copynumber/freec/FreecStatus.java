package com.hartwig.hmftools.common.copynumber.freec;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public enum FreecStatus {
    GERMLINE,
    SOMATIC,
    UNKNOWN;

    private static final double GERMLINE_MIN_RATIO = 0.75;
    private static final double GERMLINE_MAX_RATIO = 1.25;

    @NotNull
    public static FreecStatus fromNormalRatio(final double ratio) {
        if (Doubles.isZero(ratio)) {
            return UNKNOWN;
        }

        if (Doubles.greaterThan(ratio, GERMLINE_MAX_RATIO) || Doubles.lessThan(ratio, GERMLINE_MIN_RATIO)) {
            return GERMLINE;
        }

        return SOMATIC;
    }
}
