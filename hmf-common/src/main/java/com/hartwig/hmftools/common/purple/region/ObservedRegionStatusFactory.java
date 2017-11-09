package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_HET_DELETION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_HOM_DELETION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_NOISE;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.SOMATIC;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.UNKNOWN;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;

import org.jetbrains.annotations.NotNull;

class ObservedRegionStatusFactory {
    private static final double GERMLINE_HOM_DELETION_THRESHOLD = 0.2;
    private static final double GERMLINE_HET_DELETION_THRESHOLD = 0.8;
    private static final double GERMLINE_AMPLIFICATION_THRESHOLD = 1.2;
    private static final double GERMLINE_NOISE_THRESHOLD = 2.2;

    private final Gender gender;

    ObservedRegionStatusFactory(final Gender gender) {
        this.gender = gender;
    }

    ObservedRegionStatus status(@NotNull final PurpleSegment segment, final double normalRatio) {
        if (segment.svCluster()) {
            return ObservedRegionStatus.CLUSTER;
        }

        return fromNormalRatio(segment.chromosome(), normalRatio);
    }

    @NotNull
    ObservedRegionStatus fromNormalRatio(final String chromosome, final double ratio) {
        if (Doubles.isZero(ratio)) {
            return UNKNOWN;
        }
        double adjustment = chromosome.equals("X") && gender.equals(Gender.MALE) || chromosome.equals("Y") ? 2 : 1;

        if (Doubles.lessThan(ratio, GERMLINE_HOM_DELETION_THRESHOLD / adjustment)) {
            return GERMLINE_HOM_DELETION;
        }

        if (Doubles.lessThan(ratio, GERMLINE_HET_DELETION_THRESHOLD / adjustment)) {
            return GERMLINE_HET_DELETION;
        }

        if (Doubles.greaterThan(ratio, GERMLINE_NOISE_THRESHOLD / adjustment)) {
            return GERMLINE_NOISE;
        }

        if (Doubles.greaterThan(ratio, GERMLINE_AMPLIFICATION_THRESHOLD / adjustment)) {
            return GERMLINE_AMPLIFICATION;
        }

        return SOMATIC;
    }
}
