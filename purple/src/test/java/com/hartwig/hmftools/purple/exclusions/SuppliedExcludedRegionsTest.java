package com.hartwig.hmftools.purple.exclusions;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class SuppliedExcludedRegionsTest
{
    @Test
    public void testExclusion()
    {
        List<ChrBaseRegion> regions = new ArrayList<>();
        regions.add(region(_1, 100, 120));
        regions.add(region(_1, 200, 220));
        regions.add(region(_1, 1200, 1220));
        regions.add(region(_1, 13000, 15000));
        regions.add(region(_2, 11000, 12020));

        SuppliedExcludedRegions excludedRegions = new SuppliedExcludedRegions(regions);
        Map<Chromosome, List<CobaltRatio>> ratios = new HashMap<>();
        List<CobaltRatio> ratios1 = new ArrayList<>();
        ratios1.add(ratio(_1, 1, 0.11));
        ratios1.add(ratio(_1, 1001, 0.12));
        ratios1.add(ratio(_1, 2001, 0.13));
        ratios1.add(ratio(_1, 3001, 0.12));
        ratios1.add(ratio(_1, 9001, 0.14));
        ratios1.add(ratio(_1, 10001, 0.12));
        ratios1.add(ratio(_1, 11001, 0.14));
        ratios1.add(ratio(_1, 12001, 0.15));
        ratios1.add(ratio(_1, 13001, 0.13));
        ratios1.add(ratio(_1, 14001, 0.12));
        ratios1.add(ratio(_1, 15001, 0.10));
        ratios1.add(ratio(_1, 16001, 0.13));
        ratios1.add(ratio(_1, 17001, 0.15));

        List<CobaltRatio> ratios2 = new ArrayList<>();
        ratios2.add(ratio(_2, 10001, 0.12));
        ratios2.add(ratio(_2, 11001, 0.14));
        ratios2.add(ratio(_2, 12001, 0.15));

        List<CobaltRatio> ratios4 = new ArrayList<>();
        ratios4.add(ratio(_4, 10001, 0.12));
        ratios4.add(ratio(_4, 11001, 0.14));

        List<CobaltRatio> ratiosX = new ArrayList<>();
        ratiosX.add(ratio(_X, 10001, 0.12));
        ratiosX.add(ratio(_X, 11001, 0.14));
        ratiosX.add(ratio(_X, 12001, 0.15));

        ratios.put(_1, ratios1);
        ratios.put(_2, ratios2);
        ratios.put(_3, new ArrayList<>());
        ratios.put(_4, ratios4);
        ratios.put(_X, ratiosX);

        Map<Chromosome, List<CobaltRatio>> filtered = excludedRegions.maskIntersectingRatios(ratios);
        List<CobaltRatio> filtered1 = filtered.get(_1);
        assertEquals(ratios1.size(), filtered1.size());
        assertEquals(ratios1.get(0).mask(), filtered1.get(0));
        assertEquals(ratios1.get(1).mask(), filtered1.get(1));
        assertEquals(ratios1.get(2), filtered1.get(2));
        assertEquals(ratios1.get(3), filtered1.get(3));
        assertEquals(ratios1.get(4), filtered1.get(4));
        assertEquals(ratios1.get(5), filtered1.get(5));
        assertEquals(ratios1.get(6), filtered1.get(6)); // 11001-12000, not masked
        assertEquals(ratios1.get(7).mask(), filtered1.get(7)); // 12001-13000, masked
        assertEquals(ratios1.get(8).mask(), filtered1.get(8)); // 13001-14000, masked
        assertEquals(ratios1.get(9).mask(), filtered1.get(9)); // 14001-15000, masked
        assertEquals(ratios1.get(10), filtered1.get(10));
        assertEquals(ratios1.get(11), filtered1.get(11));
        assertEquals(ratios1.get(12), filtered1.get(12));

        List<CobaltRatio> filtered2 = filtered.get(_2);
        assertEquals(ratios2.size(), filtered2.size());
        assertEquals(ratios2.get(0).mask(), filtered2.get(0));
        assertEquals(ratios2.get(1).mask(), filtered2.get(1));
        assertEquals(ratios2.get(2).mask(), filtered2.get(2));

        assertEquals(0, filtered.get(_3).size());
        assertEquals(ratios4, filtered.get(_4));
        assertEquals(ratiosX, filtered.get(_X));
    }

    private ChrBaseRegion region(HumanChromosome chromosome, int start, int end)
    {
        return new ChrBaseRegion(chromosome.shortName(), start, end);
    }

    private CobaltRatio ratio(HumanChromosome chromosome, int position, double v)
    {
        return new CobaltRatio(chromosome.shortName(), position, v, v, v, v, v, v, v);
    }
}
