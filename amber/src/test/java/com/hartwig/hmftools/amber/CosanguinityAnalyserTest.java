package com.hartwig.hmftools.amber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.ArrayList;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class CosanguinityAnalyserTest
{
    @Test
    public void testFindUniparentalDisomy()
    {
        var rohs = new ArrayList<RegionOfHomozygosity>();
        rohs.add(new RegionOfHomozygosity(HumanChromosome._1, 1, 5_000_000, 0, 0, 0));
        rohs.add(new RegionOfHomozygosity(HumanChromosome._1, 1_000_000, 6_000_000, 0, 0, 0));

        Chromosome uniparentalDisomy = ConsanguinityAnalyser.findUniparentalDisomy(rohs);

        assertEquals(uniparentalDisomy, HumanChromosome._1);
    }

    // test the case where two chromosomes contain long ROHs
    @Test
    public void testFindUniparentalDisomy2()
    {
        var rohs = new ArrayList<RegionOfHomozygosity>();
        rohs.add(new RegionOfHomozygosity(HumanChromosome._1, 1, 5_000_000, 0, 0, 0));
        rohs.add(new RegionOfHomozygosity(HumanChromosome._1, 1_000_000, 6_000_000, 0, 0, 0));
        rohs.add(new RegionOfHomozygosity(HumanChromosome._2, 1, 5_000_000, 0, 0, 0));

        Chromosome uniparentalDisomy = ConsanguinityAnalyser.findUniparentalDisomy(rohs);

        assertNull(uniparentalDisomy);
    }
}
