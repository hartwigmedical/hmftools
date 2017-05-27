package com.hartwig.hmftools.common.copynumber.freec;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;

import org.jetbrains.annotations.NotNull;

public enum FreecRatioRegions {
    ;

    @NotNull
    public static List<GenomeRegion> createRegionsFromRatios(@NotNull final List<FreecRatio> ratios) {
        final List<GenomeRegion> result = Lists.newArrayList();

        FreecRatio start = null;
        FreecRatio end = null;
        for (final FreecRatio ratio : ratios) {
            if (start == null) {
                start = ratio;
                end = ratio;
            } else if (isSameRegion(start, ratio)) {
                end = ratio;
            } else {
                result.add(create(start, end.position() - 1));
                start = ratio;
                end = ratio;
            }
        }

        if (start != null) {
            result.add(create(start, end.position() + 999));
        }

        return result;
    }

    @NotNull
    private static GenomeRegion create(@NotNull final FreecRatio startOfCurrentRegion,long end) {


        return ImmutableBEDGenomeRegion.builder().chromosome(startOfCurrentRegion.chromosome()).start(startOfCurrentRegion.position()).end(
                end).build();
    }

    private static boolean isSameRegion(@NotNull final FreecRatio first, @NotNull final FreecRatio second) {
        return first.chromosome().equals(second.chromosome()) && Doubles.equal(first.medianRatio(),
                second.medianRatio()) && Doubles.equal(first.estimatedBAF(), second.estimatedBAF());
    }
}
