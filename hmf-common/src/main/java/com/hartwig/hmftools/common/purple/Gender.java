package com.hartwig.hmftools.common.purple;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public enum Gender {
    MALE,
    FEMALE,
    @Deprecated MALE_KLINEFELTER;

    private static final int MIN_BAF_COUNT = 1000;

    @NotNull
    public static Gender fromAmber(@NotNull final Multimap<Chromosome, AmberBAF> bafs) {
        return bafs.get(HumanChromosome._X).stream().filter(x -> x.position() > 2_699_520 && x.position() < 155_260_560).count()
                > MIN_BAF_COUNT ? FEMALE : MALE;
    }

}
