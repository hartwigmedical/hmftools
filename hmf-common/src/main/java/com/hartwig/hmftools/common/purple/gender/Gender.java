package com.hartwig.hmftools.common.purple.gender;

import java.util.Collection;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;

public enum Gender {
    MALE,
    FEMALE;

    private static final int MIN_BAF_COUNT = 1000;

    public static Gender fromObservedRegions(Collection<ObservedRegion> regions) {
        return regions.stream().filter(x -> x.chromosome().equals("X"))
                .filter(x -> x.end() > 2_699_520 && x.start() < 155_260_560)
                .mapToInt(ObservedRegion::bafCount).sum() > MIN_BAF_COUNT
                ? FEMALE
                : MALE;
    }

    public static Gender fromCopyNumbers(Collection<PurpleCopyNumber> regions) {
        return regions.stream().filter(x -> x.chromosome().equals("X"))
                .filter(x -> x.end() > 2_699_520 && x.start() < 155_260_560)
                .mapToInt(PurpleCopyNumber::bafCount).sum() > MIN_BAF_COUNT
                ? FEMALE
                : MALE;
    }
}
