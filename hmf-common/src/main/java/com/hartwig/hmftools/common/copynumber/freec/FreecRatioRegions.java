package com.hartwig.hmftools.common.copynumber.freec;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;

import org.jetbrains.annotations.NotNull;

public enum FreecRatioRegions {
    ;

    private static final int SEGMENT_SIZE = 1000;

    @NotNull
    public static List<GenomeRegion> createRegionsFromRatios(@NotNull final List<FreecRatio> ratios) {
        final List<GenomeRegion> result = Lists.newArrayList();

        FreecRatio start = null;
        long endPosition = 0;
        for (final FreecRatio ratio : ratios) {
            if (start == null) {
                start = ratio;
            } else if (isNewRegion(start, ratio)) {
                result.add(create(start, endPosition));
                start = ratio;
            }

            endPosition = ratio.position() + SEGMENT_SIZE - 1;
        }

        if (start != null) {
            result.add(create(start, endPosition));
        }

        return result;
    }

    @NotNull
    private static GenomeRegion create(@NotNull final FreecRatio startOfCurrentRegion,long end) {
        return ImmutableBEDGenomeRegion.builder().chromosome(startOfCurrentRegion.chromosome())
                .start(startOfCurrentRegion.position()).end(end).build();
    }

    private static boolean isNewRegion(@NotNull final FreecRatio first, @NotNull final FreecRatio second) {
        return isDiscontiguous(first, second) || !Doubles.equal(first.medianRatio(),
                second.medianRatio()) || !Doubles.equal(first.estimatedBAF(), second.estimatedBAF()) ;
    }

    private static boolean isDiscontiguous(@NotNull final FreecRatio first, @NotNull final FreecRatio second) {
        return !first.chromosome().equals(second.chromosome()) || second.position() - first.position() != 1000;
    }
}
