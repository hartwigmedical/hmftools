package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public enum ObservedRegionStatus {
    GERMLINE,
    SOMATIC,
    UNKNOWN;

    private static final double GERMLINE_MIN_RATIO = 0.8;
    private static final double GERMLINE_MAX_RATIO = 1.2;

    @NotNull
    public static ObservedRegionStatus fromNormalRatio(final Gender gender, final String chromosome, final double ratio) {
        if (Doubles.isZero(ratio)) {
            return UNKNOWN;
        }
        double adjustment = chromosome.equals("X") && gender.equals(Gender.MALE) || chromosome.equals("Y") ? 0.5 : 0;
        if (Doubles.greaterThan(ratio, GERMLINE_MAX_RATIO - adjustment) || Doubles.lessThan(ratio, GERMLINE_MIN_RATIO - adjustment)) {
            return GERMLINE;
        }

        return SOMATIC;
    }
}
