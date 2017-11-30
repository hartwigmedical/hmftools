package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_HET_DELETION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_HOM_DELETION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_NOISE;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.UNKNOWN;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;

import org.jetbrains.annotations.NotNull;

class ObservedRegionStatusFactory {
    private static final double GERMLINE_HOM_DELETION_THRESHOLD = 0.1;
    private static final double GERMLINE_HET_DELETION_THRESHOLD = 0.8;
    private static final double GERMLINE_AMPLIFICATION_THRESHOLD = 1.2;
    private static final double GERMLINE_NOISE_THRESHOLD = 2.2;

    private final Gender gender;

    ObservedRegionStatusFactory(final Gender gender) {
        this.gender = gender;
    }

    ObservedRegionStatus status(@NotNull final PurpleSegment segment, final double normalRatio, final double tumorRatio) {
        if (segment.svCluster()) {
            return ObservedRegionStatus.CLUSTER;
        }

        return fromRatio(segment.chromosome(), normalRatio, tumorRatio);
    }

    @NotNull
    ObservedRegionStatus fromRatio(final String chromosome, final double normalRatio, final double tumorRatio) {
        if (Doubles.isZero(normalRatio)) {
            return UNKNOWN;
        }
        double adjustment = chromosome.equals("X") && gender.equals(Gender.MALE) || chromosome.equals("Y") ? 2 : 1;

        double adjustedHomDeletionThreshold = GERMLINE_HOM_DELETION_THRESHOLD / adjustment;
        if (Doubles.lessThan(normalRatio, adjustedHomDeletionThreshold) && Doubles.lessThan(tumorRatio, adjustedHomDeletionThreshold)) {
            return GERMLINE_HOM_DELETION;
        }

        if (Doubles.lessThan(normalRatio, GERMLINE_HET_DELETION_THRESHOLD / adjustment)) {
            return GERMLINE_HET_DELETION;
        }

        if (Doubles.greaterThan(normalRatio, GERMLINE_NOISE_THRESHOLD / adjustment)) {
            return GERMLINE_NOISE;
        }

        if (Doubles.greaterThan(normalRatio, GERMLINE_AMPLIFICATION_THRESHOLD / adjustment)) {
            return GERMLINE_AMPLIFICATION;
        }

        return DIPLOID;
    }
}
