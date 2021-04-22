package com.hartwig.hmftools.cobalt.ratio;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

class DiploidRatioSupplier {

    private static final long ROLLING_MEDIAN_MAX_DISTANCE = 5_000;
    private static final long ROLLING_MEDIAN_MIN_COVERAGE = 1_000;

    private final ListMultimap<Chromosome, ReadRatio> result = ArrayListMultimap.create();

    DiploidRatioSupplier(@NotNull final CobaltChromosomes chromosomes, @NotNull final ListMultimap<Chromosome, ReadRatio> normalRatios) {

        for (CobaltChromosome cobaltChromosome : chromosomes.chromosomes()) {
            if (HumanChromosome.contains(cobaltChromosome.contig())) {
                Chromosome chromosome = HumanChromosome.fromString(cobaltChromosome.contig());
                final List<ReadRatio> ratios = normalRatios.get(chromosome);
                final List<ReadRatio> adjustedRatios;
                if (chromosome.equals(HumanChromosome._Y)) {
                    adjustedRatios = ratios;
                } else {
                    double expectedRatio = cobaltChromosome.actualRatio();
                    adjustedRatios = new DiploidRatioNormalization(expectedRatio,
                            ROLLING_MEDIAN_MAX_DISTANCE,
                            ROLLING_MEDIAN_MIN_COVERAGE,
                            ratios).get();
                }
                result.replaceValues(chromosome, adjustedRatios);
            }
        }
    }

    @NotNull
    ListMultimap<Chromosome, ReadRatio> result() {
        return result;
    }
}
