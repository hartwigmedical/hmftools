package com.hartwig.hmftools.amber.purity;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RegionsFilter
{
    private final Map<String, List<ChrBaseRegion>> RegionsToExclude = new HashMap<>();

    public RegionsFilter(List<ChrBaseRegion> regionsToExclude)
    {
        for(ChrBaseRegion region : regionsToExclude)
        {
            RegionsToExclude.computeIfAbsent(region.chromosome(), k -> new ArrayList<>()).add(region);
        }
        for(List<ChrBaseRegion> regions : RegionsToExclude.values())
        {
            Collections.sort(regions);
        }
    }

    public <T extends GenomePosition> List<T> filter(List<T> raw)
    {
        List<T> result = new ArrayList<>();
        String currentChromosome = null;
        List<ChrBaseRegion> currentRegions = null;
        int regionIndex = 0;

        for(T position : raw)
        {
            String chromosome = position.chromosome();
            if(!chromosome.equals(currentChromosome))
            {
                currentChromosome = chromosome;
                currentRegions = RegionsToExclude.get(chromosome);
                regionIndex = 0;
            }

            if(currentRegions == null)
            {
                result.add(position);
                continue;
            }

            while(regionIndex < currentRegions.size() && currentRegions.get(regionIndex).end() < position.position())
            {
                regionIndex++;
            }

            if(regionIndex >= currentRegions.size() || !currentRegions.get(regionIndex).containsPosition(position))
            {
                result.add(position);
            }
        }
        return result;
    }
}
