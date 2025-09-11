package com.hartwig.hmftools.cobalt.exclusions;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class SuppliedExcludedRegionsTest
{
    @Test
    public void findIntersections()
    {
        // Prepare a GC Profile map with chromosomes 1, 2, 3 and profile data for positions as follows:
        // 1: 1, 1001, ... , 19001
        // 2: 10001, 11001, ... , 19001
        // 3: 1, 1001, ... , 9001
        ListMultimap<Chromosome, GCProfile> ratios = ArrayListMultimap.create();
        int pos = 1;
        for(int i = 0; i < 10; i++)
        {
            ratios.put(_1, gcp(_1, pos));
            ratios.put(_3, gcp(_3, pos));
            pos += 1000;
        }
        for(int i = 0; i < 10; i++)
        {
            ratios.put(_1, gcp(_1, pos));
            ratios.put(_2, gcp(_2, pos));
            pos += 1000;
        }

        List<ChrBaseRegion> regions = new ArrayList<>();
        regions.add(region(_1, 100, 120));
        regions.add(region(_1, 200, 220));
        regions.add(region(_1, 1200, 1220));
        regions.add(region(_1, 13000, 15000));
        regions.add(region(_2, 1100, 1220));
        regions.add(region(_2, 11000, 12020));

        // These regions should collect:
        // 1: 1-1000, 1001-2000, 12_001-13_000, 13_001-14_000, 14_001-15_000
        // 2: 10_001-11_000, 11_001-12_000, 12_001-13_000
        SuppliedExcludedRegions excludedRegions = new SuppliedExcludedRegions(regions);
        ListMultimap<Chromosome, GCProfile> filtered = excludedRegions.findIntersections(ratios);
        assertEquals(8, filtered.size());
        List<GCProfile> filtered1 = filtered.get(_1);
        assertEquals(5, filtered1.size());
        assertEquals(gcp(_1, 1), filtered1.get(0));
        assertEquals(gcp(_1, 1_001), filtered1.get(1));
        assertEquals(gcp(_1, 12_001), filtered1.get(2));
        assertEquals(gcp(_1, 13_001), filtered1.get(3));
        assertEquals(gcp(_1, 14_001), filtered1.get(4));

        List<GCProfile> filtered2 = filtered.get(_2);
        assertEquals(3, filtered2.size());
        assertEquals(gcp(_2, 10_001), filtered2.get(0));
        assertEquals(gcp(_2, 11_001), filtered2.get(1));
        assertEquals(gcp(_2, 12_001), filtered2.get(2));
    }

    private ChrBaseRegion region(HumanChromosome chromosome, int start, int end)
    {
        return new ChrBaseRegion(chromosome.shortName(), start, end);
    }

    private GCProfile gcp(HumanChromosome chromosome, int position)
    {
        return ImmutableGCProfile.builder()
                .chromosome(chromosome.shortName())
                .start(position)
                .end(position + 1000)
                .gcContent(0.50)
                .nonNPercentage(90)
                .mappablePercentage(90)
                .build();
    }
}
