package com.hartwig.hmftools.cobalt.ratio;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.cobalt.ReferenceRatioStatistics;
import com.hartwig.hmftools.common.cobalt.ReferenceRatioStatisticsFactory;

import org.jetbrains.annotations.NotNull;

class DiploidRatioSupplier {

    private static final long ROLLING_MEDIAN_MAX_DISTANCE = 5_000;
    private static final long ROLLING_MEDIAN_MIN_COVERAGE = 1_000;

    private final ListMultimap<Chromosome, ReadRatio> result = ArrayListMultimap.create();

    DiploidRatioSupplier(@NotNull final ListMultimap<Chromosome, ReadRatio> normalRatios) {

        final ReferenceRatioStatistics stats = ReferenceRatioStatisticsFactory.fromReferenceRatio(normalRatios);

        for (Chromosome chromosome : normalRatios.keySet()) {
            final List<ReadRatio> ratios = normalRatios.get(chromosome);
            final List<ReadRatio> adjustedRatios;
            if (chromosome.equals(HumanChromosome._Y)) {
                adjustedRatios = ratios;
            } else {
                double expectedRatio = chromosome.equals(HumanChromosome._X) && !stats.containsTwoXChromosomes() ? 0.5 : 1;
                adjustedRatios = new DiploidRatioNormalization(expectedRatio,
                        ROLLING_MEDIAN_MAX_DISTANCE,
                        ROLLING_MEDIAN_MIN_COVERAGE,
                        ratios).get();
            }
            result.replaceValues(chromosome, adjustedRatios);
        }
    }

    @NotNull
    ListMultimap<Chromosome, ReadRatio> result() {
        return result;
    }
}
