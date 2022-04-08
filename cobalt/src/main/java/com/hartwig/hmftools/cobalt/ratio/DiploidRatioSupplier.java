package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import java.util.Collection;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public final class DiploidRatioSupplier
{
    public static ArrayListMultimap<Chromosome, ReadRatio> calcDiploidRatioResults(
            final Collection<Chromosome> chromosomeList, final ListMultimap<Chromosome, ReadRatio> normalRatios, final List<MedianRatio> medianRatios)
    {
        ArrayListMultimap<Chromosome, ReadRatio> results = ArrayListMultimap.create();

        for (CobaltChromosome cobaltChromosome : new CobaltChromosomes(medianRatios).chromosomes())
        {
            Chromosome chr = Chromosome.findByContig(cobaltChromosome.contig(), chromosomeList);
            if(HumanChromosome.contains(cobaltChromosome.contig()))
            {
                final List<ReadRatio> ratios = normalRatios.get(chr);
                final List<ReadRatio> adjustedRatios;
                if (HumanChromosome.fromString(cobaltChromosome.contig()).equals(HumanChromosome._Y))
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
                results.replaceValues(chr, adjustedRatios);
            }
        }

        return results;
    }
}
