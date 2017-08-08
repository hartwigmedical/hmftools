package com.hartwig.hmftools.common.purple.gender;

import java.util.Collection;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.baf.TumorBAF;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

public enum Gender {
    MALE,
    FEMALE;

    private static final int MIN_BAF_COUNT = 1000;

    public static Gender fromBAFCount(@NotNull final Multimap<String, TumorBAF> bafs) {
        return bafs.get("X").stream().filter(x -> x.position() > 2_699_520 && x.position() < 155_260_560).count() > MIN_BAF_COUNT
                ? FEMALE
                : MALE;
    }

    public static Gender fromObservedRegions(Collection<ObservedRegion> regions) {
        return regions.stream()
                .filter(x -> x.chromosome().equals("X"))
                .filter(x -> x.end() > 2_699_520 && x.start() < 155_260_560)
                .mapToInt(ObservedRegion::bafCount)
                .sum() > MIN_BAF_COUNT ? FEMALE : MALE;
    }

    public static Gender fromCopyNumbers(Collection<PurpleCopyNumber> regions) {
        return regions.stream()
                .filter(x -> x.chromosome().equals("X"))
                .filter(x -> x.end() > 2_699_520 && x.start() < 155_260_560)
                .mapToInt(PurpleCopyNumber::bafCount)
                .sum() > MIN_BAF_COUNT ? FEMALE : MALE;
    }
}
