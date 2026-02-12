package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Set;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;

class GermlineGainDeletionUtil
{
    public static double getSomaticMinCopyNumber(final List<GermlineAmpDel> gainDelsForGene)
    {
        double minCopyNumberFromDeletions = getMinimumTumorCopyNumberFromAmpDels(gainDelsForGene);
        return Math.max(0, minCopyNumberFromDeletions);
    }

    public static String getChromosome(final List<GermlineAmpDel> deletions)
    {
        if(deletions.isEmpty())
        {
            throw new IllegalArgumentException("No germline deletions provided to get chromosome from.");
        }

        Set<String> chromosomes = deletions.stream().map(d -> d.Chromosome).collect(Collectors.toSet());
        if(chromosomes.size() == 1)
        {
            return chromosomes.iterator().next();
        }
        else
        {
            throw new IllegalStateException("Germline deletions disagree on chromosome.");
        }
    }

    public static String getChromosomeBand(final List<GermlineAmpDel> deletions)
    {
        if(deletions.isEmpty())
        {
            throw new IllegalArgumentException("No germline deletions provided to get chromosome band from.");
        }

        Set<String> chromosomeBands = deletions.stream().map(d -> d.ChromosomeBand).collect(Collectors.toSet());
        if(chromosomeBands.size() == 1)
        {
            return chromosomeBands.iterator().next();
        }
        else
        {
            throw new IllegalStateException("Germline deletions disagree on chromosome band.");
        }
    }

    public static GeneCopyNumber findGeneCopyNumberForGene(final String geneNameToFind, final List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        return allSomaticGeneCopyNumbers.stream()
                .filter(g -> g.geneName().equals(geneNameToFind))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("Could not find gene copy number for gene with name: " + geneNameToFind));
    }

    private static double getMinimumTumorCopyNumberFromAmpDels(final List<GermlineAmpDel> gainDelsForGene)
    {
        return getCopyNumberFromAmpDels(gainDelsForGene,
                germlineAmpDel -> germlineAmpDel.TumorCopyNumber,
                Double::min);
    }

    private static double getCopyNumberFromAmpDels(
            final List<GermlineAmpDel> ampDelsForGene,
            Function<GermlineAmpDel, Double> copyNumberGetter,
            BinaryOperator<Double> minMaxOperator)
    {
        if(ampDelsForGene.isEmpty())
        {
            throw new IllegalArgumentException("Cannot determine copy number for empty list of germline amp dels");
        }

        return Math.max(0, ampDelsForGene.stream().map(copyNumberGetter).reduce(minMaxOperator).orElseThrow());
    }
}
