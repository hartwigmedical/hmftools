package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;

import org.jetbrains.annotations.NotNull;

class ObservedRegionStatusFactory {
    private final Gender gender;

    ObservedRegionStatusFactory(final Gender gender) {
        this.gender = gender;
    }

    ObservedRegionStatus status(@NotNull final PurpleSegment segment, final double normalRatio) {
        switch (segment.status()) {
            case CENTROMERE:
                return ObservedRegionStatus.CENTROMERE;
            case CLUSTER:
                return ObservedRegionStatus.CLUSTER;
        }
        return ObservedRegionStatus.fromNormalRatio(gender, segment.chromosome(), normalRatio);
    }

}
