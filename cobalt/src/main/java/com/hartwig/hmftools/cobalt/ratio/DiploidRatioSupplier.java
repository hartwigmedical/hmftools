package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public final class DiploidRatioSupplier
{
    public static ListMultimap<Chromosome, ReadRatio> calcDiploidRatioResults(
            final CobaltChromosomes chromosomes, final ListMultimap<Chromosome, ReadRatio> normalRatios)
    {
        ListMultimap<Chromosome, ReadRatio> results = ArrayListMultimap.create();

        for(CobaltChromosome cobaltChromosome : chromosomes.chromosomes())
        {
            if(HumanChromosome.contains(cobaltChromosome.contig()))
            {
                Chromosome chromosome = HumanChromosome.fromString(cobaltChromosome.contig());
                final List<ReadRatio> ratios = normalRatios.get(chromosome);
                final List<ReadRatio> adjustedRatios;
                if(chromosome.equals(HumanChromosome._Y))
                {
                    adjustedRatios = ratios;
                }
                else
                {
                    double expectedRatio = cobaltChromosome.actualRatio();
                    adjustedRatios = new DiploidRatioNormalization(expectedRatio,
                            ROLLING_MEDIAN_MAX_DISTANCE,
                            ROLLING_MEDIAN_MIN_COVERAGE,
                            ratios).get();
                }

                results.replaceValues(chromosome, adjustedRatios);
            }
        }

        return results;
    }
}
