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
        FreecRatio end = null;
        for (final FreecRatio ratio : ratios) {
            if (start == null) {
                start = ratio;
            } else if (isNewRegion(start, ratio)) {
                result.add(create(start, end));
                start = ratio;
            }
            end = ratio;
        }

        if (start != null) {
            result.add(create(start, end));
        }

        return result;
    }

    @NotNull
    private static GenomeRegion create(@NotNull final FreecRatio start, FreecRatio end) {
        return ImmutableBEDGenomeRegion.builder().chromosome(start.chromosome())
                .start(start.position()).end(end.position() + SEGMENT_SIZE - 1).build();
    }

    private static boolean isNewRegion(@NotNull final FreecRatio first, @NotNull final FreecRatio second) {
        return !first.chromosome().equals(second.chromosome()) || !Doubles.equal(first.medianRatio(),
                second.medianRatio()) || !Doubles.equal(first.estimatedBAF(), second.estimatedBAF());
    }

}
