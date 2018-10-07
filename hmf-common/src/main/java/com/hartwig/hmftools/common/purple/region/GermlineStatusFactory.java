package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.purple.region.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.NOISE;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.UNKNOWN;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;

import org.jetbrains.annotations.NotNull;

class GermlineStatusFactory {
    private static final double GERMLINE_HOM_DELETION_THRESHOLD = 0.1;
    private static final double GERMLINE_HET_DELETION_THRESHOLD = 0.85;
    private static final double GERMLINE_AMPLIFICATION_THRESHOLD = 1.15;
    private static final double GERMLINE_NOISE_THRESHOLD = 2.2;

    private final Gender gender;

    GermlineStatusFactory(final Gender gender) {
        this.gender = gender;
    }

    GermlineStatus status(@NotNull final PurpleSegment segment, final double normalRatio, final double tumorRatio) {
        return fromRatio(segment.chromosome(), normalRatio, tumorRatio);
    }

    @NotNull
    GermlineStatus fromRatio(final String chromosome, final double normalRatio, final double tumorRatio) {
        if (Doubles.isZero(normalRatio)) {
            return UNKNOWN;
        }
        double adjustment = chromosome.equals("X") && gender.equals(Gender.MALE) || chromosome.equals("Y") ? 2 : 1;

        double adjustedHomDeletionThreshold = GERMLINE_HOM_DELETION_THRESHOLD / adjustment;
        if (Doubles.lessThan(normalRatio, adjustedHomDeletionThreshold) && Doubles.lessThan(tumorRatio, adjustedHomDeletionThreshold)) {
            return HOM_DELETION;
        }

        if (Doubles.lessThan(normalRatio, GERMLINE_HET_DELETION_THRESHOLD / adjustment)) {
            return HET_DELETION;
        }

        if (Doubles.greaterThan(normalRatio, GERMLINE_NOISE_THRESHOLD / adjustment)) {
            return NOISE;
        }

        if (Doubles.greaterThan(normalRatio, GERMLINE_AMPLIFICATION_THRESHOLD / adjustment)) {
            return AMPLIFICATION;
        }

        return DIPLOID;
    }
}
