package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.FittingRegion;

import org.junit.Before;
import org.junit.Test;

public class AmberPointsProviderTest
{
    private ListMultimap<Chromosome, AmberBAF> AmberData;
    private AmberPointsProvider provider;

    @Before
    public void setup()
    {
        AmberData = ArrayListMultimap.create();
        provider = null;
    }

    @Test
    public void noPointsInRegionTest()
    {
        int position = 1_000_000;
        for(int i = 0; i < 10; i++)
        {
            addDataPoint(_2, position);
            position += 1000;
        }
        createProvider();
        assertTrue(provider.nextBlock(fr(_2, 500_000, 501_000)).isEmpty()); // before the points
        assertTrue(provider.nextBlock(fr(_2, 1_500_000, 1_501_000)).isEmpty()); // after the points
        assertTrue(provider.nextBlock(fr(_2, position - 900, position - 800)).isEmpty()); // between points
    }

    @Test
    public void regionContainsPointsTest()
    {
        int start = 1_000_000;
        int position = start;
        final int step = 1000;
        for(int i = 0; i < 10; i++)
        {
            addDataPoint(_2, position);
            position += step;
        }
        createProvider();

        FittingRegion region = fr(_2, start - 10, start + 10);
        assertEquals(1, provider.nextBlock(region).size());

        region = fr(_2, start, start + 1);
        assertEquals(1, provider.nextBlock(region).size());

        region = fr(_2, start, start + step * 5 - 1);
        assertEquals(5, provider.nextBlock(region).size());
    }

    @Test
    public void successiveRegionsAndChromosomesTest()
    {
        int start = 1_000_000;
        int position = start;
        final int step = 1000;
        for(int i = 0; i < 100; i++)
        {
            addDataPoint(_2, position);
            addDataPoint(_3, position);
            addDataPoint(_4, position);
            position += step;
        }
        createProvider();

        FittingRegion region = fr(_2, start - 1010, start - 10);
        assertEquals(0, provider.nextBlock(region).size());

        region = fr(_2, start + 10, start + 10 + step * 5);
        assertEquals(5, provider.nextBlock(region).size());

        region = fr(_2, start + 10 + step * 5 + 10, start + 10 + step * 5 + 110); // small region between points
        assertEquals(0, provider.nextBlock(region).size());

        region = fr(_2, start + 10 + step * 10, start + 10 + step * 20);
        assertEquals(10, provider.nextBlock(region).size());

        region = fr(_2, start + 10 + step * 30, start + 10 + step * 60);
        assertEquals(30, provider.nextBlock(region).size());

        region = fr(_2, start + 10 + step * 60 + 10, start + 10 + step * 60 + 110); // small region between points
        assertEquals(0, provider.nextBlock(region).size());

        region = fr(_2, start + 10 + step * 70, start + 10 + step * 75);
        assertEquals(5, provider.nextBlock(region).size());

        region = fr(_2, start + 10 + step * 100, start + 10 + step * 105);
        assertEquals(0, provider.nextBlock(region).size());

        region = fr(_3, start + 10, start + 10 + step * 5);
        assertEquals(5, provider.nextBlock(region).size());

        region = fr(_3, start + 10 + step * 10, start + 10 + step * 20);
        assertEquals(10, provider.nextBlock(region).size());

        region = fr(_4, start + 10, start + 10 + step * 5);
        assertEquals(5, provider.nextBlock(region).size());

        region = fr(_4, start + 10 + step * 10, start + 10 + step * 20);
        assertEquals(10, provider.nextBlock(region).size());
    }

    private void createProvider()
    {
        provider = new AmberPointsProvider(AmberData);
    }

    private void addDataPoint(HumanChromosome chromosome, int position)
    {
        AmberData.put(chromosome, ab(chromosome, position));
    }

    private AmberBAF ab(HumanChromosome chromosome, int position)
    {
        return new AmberBAF(V38.versionedChromosome(chromosome), position, 0.5, 100, -1, 0);
    }

    private FittingRegion fr(HumanChromosome chromosome, int start, int end)
    {
        return new FR(chromosome, start, end, GermlineStatus.DIPLOID, 10, 0.5, 0.5, 0.5);
    }
}
