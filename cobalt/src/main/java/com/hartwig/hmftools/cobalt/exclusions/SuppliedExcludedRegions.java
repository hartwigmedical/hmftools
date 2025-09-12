package com.hartwig.hmftools.cobalt.exclusions;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toList;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class SuppliedExcludedRegions
{
    private final Map<Chromosome, List<ChrBaseRegion>> regionsByChromosome = new HashMap<>();

    public SuppliedExcludedRegions(final List<ChrBaseRegion> regions)
    {
        Map<String, List<ChrBaseRegion>> grouped = regions.stream().collect(groupingBy(ChrBaseRegion::chromosome));
        grouped.forEach((key, value) ->
                regionsByChromosome.put(HumanChromosome.fromString(key), value.stream().sorted().collect(toList())));
    }

    public ListMultimap<Chromosome, GCProfile> findIntersections(ListMultimap<Chromosome, GCProfile> gcData)
    {
        ListMultimap<Chromosome, GCProfile> result = ArrayListMultimap.create();
        gcData.keySet().forEach(chromosome ->
        {
            List<ChrBaseRegion> filterRegions = regionsByChromosome.get(chromosome);
            if(filterRegions != null)
            {
                List<GCProfile> gcProfilesForChromosome = gcData.get(chromosome); // sorted already

                int indexOfLastAffectedProfile = -1;
                for(ChrBaseRegion region : filterRegions)
                {
                    boolean keepSearching = true;
                    boolean found = false;
                    for(int j = indexOfLastAffectedProfile + 1; j < gcProfilesForChromosome.size() && keepSearching; j++)
                    {
                        GCProfile profile = gcProfilesForChromosome.get(j);
                        if(region.overlaps(profile.chrBaseRegion()))
                        {
                            result.put(chromosome, profile);
                            indexOfLastAffectedProfile = j;
                            found = true;
                        }
                        else
                        {
                            keepSearching = !found;
                        }
                    }
                }
            }
        });

        return result;
    }
}
