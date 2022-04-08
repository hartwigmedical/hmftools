package com.hartwig.hmftools.cobalt;

import java.util.Collection;
import java.util.function.Function;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class CobaltUtils
{
    public static <T> Multimap<com.hartwig.hmftools.common.genome.chromosome.Chromosome, T> toCommonChromosomeMap(
            Multimap<Chromosome, T> input)
    {
        Multimap<com.hartwig.hmftools.common.genome.chromosome.Chromosome, T> output = ArrayListMultimap.create();
        for (Chromosome c : input.keySet())
        {
            if (HumanChromosome.contains(c.contig))
                output.putAll(HumanChromosome.fromString(c.contig), input.get(c));
        }
        return output;
    }

    public static <T> Multimap<Chromosome, T> toChromosomeMultiMap(
            Collection<T> genomePositions, Collection<Chromosome> chromosomes,
            Function<T, String> chromosomeGetter)
    {
        Multimap<Chromosome, T> result = ArrayListMultimap.create();
        for (T gp : genomePositions)
        {
            Chromosome c = Chromosome.findByContig(chromosomeGetter.apply(gp), chromosomes);
            result.put(c, gp);
        }
        return result;
    }
}
