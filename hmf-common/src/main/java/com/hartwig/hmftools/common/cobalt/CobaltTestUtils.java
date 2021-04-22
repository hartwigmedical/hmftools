package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.MIN_Y_COUNT;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class CobaltTestUtils
{
    @NotNull
    public static CobaltChromosomes female() {
        List<MedianRatio> ratios = Lists.newArrayList();
        for (int i = 0; i < 22; i++) {
            ratios.add(create(String.valueOf(i), 1, 1));
        }

        ratios.add(create("X", 1, 1));
        return new CobaltChromosomes(ratios);
    }

    @NotNull
    public static CobaltChromosomes male() {
        List<MedianRatio> ratios = Lists.newArrayList();
        for (int i = 0; i < 22; i++) {
            ratios.add(create(String.valueOf(i), 1, 1));
        }

        ratios.add(create("X", 0.5, 1));
        ratios.add(create("Y", 0.5, MIN_Y_COUNT));
        return new CobaltChromosomes(ratios);
    }

    @NotNull
    public static MedianRatio create(@NotNull String contig, double ratio) {
        return create(contig, ratio, MIN_Y_COUNT);
    }

    @NotNull
    public static MedianRatio create(@NotNull String contig, double ratio, int count) {
        return ImmutableMedianRatio.builder().count(count).chromosome(contig).medianRatio(ratio).build();
    }

}
