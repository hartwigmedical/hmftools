package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;

import org.junit.Assert;
import org.junit.Test;

public class ChromosomeRegionCountsTest extends ToolsTestBase
{
    @Test
    public void constructorTest()
    {
        ChromosomeRegionCounts crCounts = new ChromosomeRegionCounts(_1, new ArrayList<>());
        Assert.assertTrue(crCounts.counts().isEmpty());

        RegionGeneEvents rge1 = rge("G", _1, 1000, 2000, AMPLIFICATION);
        RegionGeneEvents rge2 = rge("G", _1, 2000, 2500, AMPLIFICATION);
        RegionGeneEvents rge3 = rge("G", _1, 2000, 3000, HET_DELETION);
        crCounts = new ChromosomeRegionCounts(_1, List.of(rge1, rge2, rge3));
        final SortedMap<RegionAmpDel, Integer> counts = crCounts.counts();
        assertEquals(3, counts.size());
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge1)));
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge2)));
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge3)));
    }

    @Test
    public void addCountsToEmptyTest()
    {
        ChromosomeRegionCounts crCounts = new ChromosomeRegionCounts(_1, new ArrayList<>());
        RegionGeneEvents rge1 = rge("G", _1, 1000, 2000, AMPLIFICATION);
        RegionGeneEvents rge2 = rge("G", _1, 2000, 2500, AMPLIFICATION);
        ChromosomeRegionCounts crCounts1 = new ChromosomeRegionCounts(_1, List.of(rge1, rge2));
        crCounts.addAll(crCounts1);

        final SortedMap<RegionAmpDel, Integer> counts = crCounts.counts();
        assertEquals(2, counts.size());
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge1)));
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge2)));
    }

    @Test
    public void addAllTest()
    {
        RegionGeneEvents rge1 = rge("G", _1, 1000, 2000, AMPLIFICATION);
        RegionGeneEvents rge2 = rge("C", _1, 2000, 2500, AMPLIFICATION);
        RegionGeneEvents rge3 = rge("G", _1, 2000, 3000, HET_DELETION);
        ChromosomeRegionCounts crCounts1 = new ChromosomeRegionCounts(_1, List.of(rge1, rge2, rge3));

        RegionGeneEvents rge4 = rge("A", _1, 2100, 2500, AMPLIFICATION);
        RegionGeneEvents rge5 = rge("B", _1, 2100, 9500, AMPLIFICATION);
        ChromosomeRegionCounts crCounts2 = new ChromosomeRegionCounts(_1, List.of(rge2, rge3, rge4, rge5));

        crCounts1.addAll(crCounts2);

        final SortedMap<RegionAmpDel, Integer> counts = crCounts1.counts();
        assertEquals(5, counts.size());
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge1)));
        assertEquals(2, (int) counts.get(new RegionAmpDel(rge2)));
        assertEquals(2, (int) counts.get(new RegionAmpDel(rge3)));
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge4)));
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge5)));
    }

    @Test
    public void regionsMustAllBeOnTheSameChromosome()
    {
        final RegionGeneEvents rge1 = rge("G", _1, 1000, 2000, AMPLIFICATION);
        final RegionGeneEvents rge2 = rge("G", _2, 1000, 2000, AMPLIFICATION);
        assertThrows(IllegalArgumentException.class, () -> new ChromosomeRegionCounts(_1, List.of(rge1, rge2)));
    }

    @Test
    public void twoEventsThatBecomeSameKindOfAmpDelTest()
    {
        ChromosomeRegionCounts crCounts = new ChromosomeRegionCounts(_1, new ArrayList<>());
        Assert.assertTrue(crCounts.counts().isEmpty());

        RegionGeneEvents rge1 = rge(_1, 1000, 2000, AMPLIFICATION);
        RegionGeneEvents rge2 = rge(_1, 1000, 2000, HOM_DELETION);
        RegionGeneEvents rge3 = rge(_1, 1000, 2000, HET_DELETION);
        crCounts = new ChromosomeRegionCounts(_1, List.of(rge1, rge2, rge3));
        final SortedMap<RegionAmpDel, Integer> counts = crCounts.counts();
        assertEquals(2, counts.size());
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge1)));
        assertEquals(2, (int) counts.get(new RegionAmpDel(rge2)));
    }

    @Test
    public void addAllWithEventsThatMergeToSameAmpDelTest()
    {
        RegionGeneEvents rge1 = rge(_1, 1000, 2000, AMPLIFICATION);
        RegionGeneEvents rge2 = rge(_1, 2000, 3000, HOM_DELETION);
        RegionGeneEvents rge3 = rge(_1, 2000, 3000, HET_DELETION);
        ChromosomeRegionCounts crCounts1 = new ChromosomeRegionCounts(_1, List.of(rge1, rge2, rge3));

        RegionGeneEvents rge4 = rge(_1, 2000, 3000, AMPLIFICATION);
        RegionGeneEvents rge5 = rge(_1, 2000, 3000, HOM_DELETION);
        ChromosomeRegionCounts crCounts2 = new ChromosomeRegionCounts(_1, List.of(rge4, rge5));

        crCounts1.addAll(crCounts2);

        final SortedMap<RegionAmpDel, Integer> counts = crCounts1.counts();
        assertEquals(3, counts.size());
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge1)));
        assertEquals(3, (int) counts.get(new RegionAmpDel(rge2)));
        assertEquals(1, (int) counts.get(new RegionAmpDel(rge4)));
    }
}
