package com.hartwig.hmftools.purple.exclusions;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class SuppliedExcludedRegions
{
    private final Map<Chromosome, List<ChrBaseRegion>> regionsByChromosome;

    public SuppliedExcludedRegions(final List<ChrBaseRegion> regions)
    {
        Map<String, List<ChrBaseRegion>> grouped = regions.stream().collect(groupingBy(ChrBaseRegion::chromosome));
        regionsByChromosome = grouped.entrySet().stream().collect(toMap(e -> HumanChromosome.fromString(e.getKey()), Map.Entry::getValue));
    }

    public Map<Chromosome, List<CobaltRatio>> maskIntersectingRatios(Map<Chromosome, List<CobaltRatio>> chromosomesToRatios)
    {
        Map<Chromosome, List<CobaltRatio>> result = new HashMap<>();
        chromosomesToRatios.forEach((chromosome, cobaltRatios) ->
        {
            List<CobaltRatio> originalRatios = chromosomesToRatios.get(chromosome);
            List<ChrBaseRegion> filterRegions = regionsByChromosome.get(chromosome);
            if(filterRegions == null)
            {
                result.put(chromosome, originalRatios);
            }
            else
            {
                List<CobaltRatio> newRatios = new ArrayList<>(originalRatios.size());
                originalRatios.forEach(cobaltRatio ->
                {
                    boolean noIntersections = cobaltRatio.findWindowOverlaps(filterRegions).isEmpty();
                    if(noIntersections)
                    {
                        newRatios.add(cobaltRatio);
                    }
                    else
                    {
                        newRatios.add(cobaltRatio.mask());
                    }
                });
                result.put(chromosome, newRatios);
            }
        });
        return result;
    }
}
