package com.hartwig.hmftools.common.cobalt;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.utils.Doubles;

public final class MedianRatioFactory
{
    public static List<MedianRatio> createFromReadRatio(final Multimap<Chromosome, ReadRatio> ratios)
    {
        return create(ReadRatio::ratio, ratios);
    }

    public static List<MedianRatio> create(final Multimap<Chromosome, CobaltRatio> ratios)
    {
        return create(CobaltRatio::referenceGCRatio, ratios);
    }

    public static List<MedianRatio> create(final Map<Chromosome,List<CobaltRatio>> chrRatiosMap)
    {
        List<MedianRatio> results = new ArrayList<>();

        for(Chromosome chromosome : chrRatiosMap.keySet())
        {
            List<CobaltRatio> ratios = chrRatiosMap.get(chromosome);

            String chromosomeStr = ratios.stream().findFirst().map(GenomePosition::chromosome).orElse(chromosome.toString());

            List<Double> contigRatios = ratios.stream()
                    .map(x -> x.referenceGCRatio()).filter(Doubles::positive).collect(Collectors.toList());

            int count = contigRatios.size();

            double medianRatio = count > 0 ? Doubles.median(contigRatios) : 0;

            results.add(new MedianRatio(chromosomeStr, medianRatio, count));
        }

        return results;
    }

    public static <T extends GenomePosition> List<MedianRatio> create(
            final Function<T, Double> ratioFunction, final Multimap<Chromosome, T> ratios)
    {
        List<MedianRatio> results = new ArrayList<>();

        for(Chromosome humanChromosome : HumanChromosome.values())
        {
            if(ratios.containsKey(humanChromosome))
            {
                Collection<T> chromosomeRatios = ratios.get(humanChromosome);
                if(!chromosomeRatios.isEmpty())
                {
                    String contig = chromosomeRatios.stream().findFirst().map(GenomePosition::chromosome).orElse(humanChromosome.toString());

                    List<Double> contigRatios = chromosomeRatios.stream().map(ratioFunction).filter(Doubles::positive).collect(Collectors.toList());

                    int count = contigRatios.size();

                    double medianRatio = count > 0 ? Doubles.median(contigRatios) : 0;

                    results.add(new MedianRatio(contig, medianRatio, count));
                }
            }
        }

        return results;
    }
}
