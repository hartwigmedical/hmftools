package com.hartwig.hmftools.common.copynumber.freec;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public enum FreecRatioRegions {
    ;

    @NotNull
    public static List<GenomeRegion> createRegionsFromRatios(@NotNull List<FreecRatio> ratios) {

        List<GenomeRegion> result = Lists.newArrayList();

        FreecRatio start = null;
        FreecRatio end = null;
        for (FreecRatio ratio : ratios) {

            if (start == null) {
                start = ratio;
                end = ratio;
            } else if (isSameRegion(start, ratio)) {
                end = ratio;
            } else {
                result.add(create(start, end));
                start = ratio;
                end = ratio;
            }
        }

        if (start != null) {
            result.add(create(start, end));
        }

        return result;
    }

    private static FreecCopyNumber create(@NotNull FreecRatio first, @NotNull FreecRatio second) {
        return ImmutableFreecCopyNumber.builder()
                .chromosome(first.chromosome())
                .start(first.position())
                .end(second.position() - 1)
                .value(0)
                .build();
    }

    private static boolean isSameRegion(@NotNull FreecRatio first, @NotNull FreecRatio second) {
        return first.chromosome().equals(second.chromosome()) &&
                Doubles.equal(first.medianRatio(), second.medianRatio()) &&
                Doubles.equal(first.estimatedBAF(), second.estimatedBAF());
    }
}
