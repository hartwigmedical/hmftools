package com.hartwig.hmftools.common.cobalt;

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

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class MedianRatioFactory
{
    @NotNull
    public static List<MedianRatio> createFromReadRatio(@NotNull Multimap<Chromosome, ReadRatio> ratios)
    {
        return create(ReadRatio::ratio, ratios);
    }

    @NotNull
    public static List<MedianRatio> create(Multimap<Chromosome, CobaltRatio> ratios)
    {
        return create(CobaltRatio::referenceGCRatio, ratios);
    }

    @NotNull
    public static List<MedianRatio> create(final Map<Chromosome,List<CobaltRatio>> chrRatiosMap)
    {
        final List<MedianRatio> results = Lists.newArrayList();

        for(Chromosome chromosome : chrRatiosMap.keySet())
        {
            final List<CobaltRatio> ratios = chrRatiosMap.get(chromosome);

            final String contig = chromosome.toString();

            final List<Double> contigRatios = ratios.stream()
                    .map(x -> x.referenceGCRatio()).filter(Doubles::positive).collect(Collectors.toList());

            int count = contigRatios.size();

            final double medianRatio = count > 0 ? Doubles.median(contigRatios) : 0;

            results.add(ImmutableMedianRatio.builder()
                    .chromosome(contig)
                    .medianRatio(medianRatio)
                    .count(count)
                    .build());
        }

        return results;
    }

    @NotNull
    public static <T extends GenomePosition> List<MedianRatio> create(@NotNull Function<T, Double> ratioFunction,
            @NotNull Multimap<Chromosome, T> ratios)
    {
        final List<MedianRatio> results = Lists.newArrayList();

        for(Chromosome humanChromosome : HumanChromosome.values())
        {
            if(ratios.containsKey(humanChromosome))
            {
                final Collection<T> chromosomeRatios = ratios.get(humanChromosome);
                if(!chromosomeRatios.isEmpty())
                {
                    final String contig =
                            chromosomeRatios.stream().findFirst().map(GenomePosition::chromosome).orElse(humanChromosome.toString());
                    final List<Double> contigRatios =
                            chromosomeRatios.stream().map(ratioFunction).filter(Doubles::positive).collect(Collectors.toList());
                    int count = contigRatios.size();
                    final double medianRatio = count > 0 ? Doubles.median(contigRatios) : 0;
                    results.add(ImmutableMedianRatio.builder()
                            .chromosome(contig)
                            .medianRatio(medianRatio)
                            .count(count)
                            .build());
                }
            }
        }
        return results;
    }
}
