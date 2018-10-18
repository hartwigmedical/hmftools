package com.hartwig.hmftools.common.purple.gender;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ReferenceRatioStatistics;
import com.hartwig.hmftools.common.cobalt.ReferenceRatioStatisticsFactory;

import org.jetbrains.annotations.NotNull;

public enum Gender {
    MALE,
    FEMALE,
    MALE_KLINEFELTER;

    private static final int MIN_BAF_COUNT = 1000;

    @NotNull
    public static Gender fromAmber(@NotNull final Multimap<Chromosome, AmberBAF> bafs) {
        return bafs.get(HumanChromosome._X).stream().filter(x -> x.position() > 2_699_520 && x.position() < 155_260_560).count()
                > MIN_BAF_COUNT ? FEMALE : MALE;
    }

    @NotNull
    public static Gender fromCobalt(@NotNull final Multimap<Chromosome, CobaltRatio> readRatios) {
        return fromCobalt(ReferenceRatioStatisticsFactory.fromCobalt(readRatios));
    }

    @NotNull
    @VisibleForTesting
    static Gender fromCobalt(@NotNull final ReferenceRatioStatistics stats) {
        if (stats.containsTwoXChromosomes()) {
            return stats.containsYChromosome() ? MALE_KLINEFELTER : FEMALE;
        }

        return MALE;
    }

}
