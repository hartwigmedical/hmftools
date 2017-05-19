package com.hartwig.hmftools.common.purple.region;

import static java.util.stream.Collectors.toList;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedCopyNumber;
import com.hartwig.hmftools.common.region.GenomeRegion;

public class ConsolidatedRegionFactory {
    public static List<ConsolidatedRegion> broad(List<FittedCopyNumber> copyNumbers) {
        return new BroadRegions().broad(copyNumbers);
    }

    public static List<ConsolidatedRegion> smooth(List<FittedCopyNumber> copyNumbers,
            List<ConsolidatedRegion> broadRegions) {

        final List<ConsolidatedRegion> result = Lists.newArrayList();

        final Set<String> orderedChromosomes = new LinkedHashSet<>();
        for (ConsolidatedRegion broadRegion : broadRegions) {
            orderedChromosomes.add(broadRegion.chromosome());
        }

        for (String orderedChromosome : orderedChromosomes) {

            final List<FittedCopyNumber> chromosomeCopyNumbers = copyNumbers.stream()
                    .filter(matchesChromosome(orderedChromosome))
                    .collect(toList());

            final List<ConsolidatedRegion> chromosomeBroadRegions = broadRegions.stream()
                    .filter(matchesChromosome(orderedChromosome))
                    .collect(toList());

            final List<ConsolidatedRegion> smoothRegions = new SmoothedRegions(chromosomeBroadRegions,
                    chromosomeCopyNumbers).getSmoothedRegions();

            result.addAll(MergeRegionBreaks.merge(smoothRegions));
        }

        return result;
    }

    private static <T extends GenomeRegion> Predicate<T> matchesChromosome(String chromosome) {
        return t -> t.chromosome().equals(chromosome);
    }

}
