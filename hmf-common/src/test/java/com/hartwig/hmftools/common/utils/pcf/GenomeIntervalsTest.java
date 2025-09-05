package com.hartwig.hmftools.common.utils.pcf;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class GenomeIntervalsTest extends PCFTestBase
{
    @Test
    public void empty()
    {
        assertEquals(0, new GenomeIntervals(new ArrayList<>()).regionsList().size());
    }

    @Test
    public void oneChromosome()
    {
        List<ChrBaseRegion> regions1 = new ArrayList<>();
        regions1.add(region(_1, 101, 200));
        PCFIntervals intervals1 = new PCFIntervals(_1, regions1);

        GenomeIntervals data = new GenomeIntervals(List.of(intervals1));
        assertEquals(regions1, data.regionsList());
    }

    @Test
    public void regionsList()
    {
        List<ChrBaseRegion> regions1 = new ArrayList<>();
        regions1.add(region(_1, 101, 200));
        regions1.add(region(_1, 201, 300));
        regions1.add(region(_1, 401, 500));
        PCFIntervals intervals1 = new PCFIntervals(_1, regions1);

        List<ChrBaseRegion> regions2 = new ArrayList<>();
        regions2.add(region(_2, 1101, 1200));
        regions2.add(region(_2, 1201, 1300));
        regions2.add(region(_2, 1401, 1500));
        PCFIntervals intervals2 = new PCFIntervals(_2, regions2);

        List<ChrBaseRegion> regions3 = new ArrayList<>();
        PCFIntervals intervals3 = new PCFIntervals(_3, regions3);

        GenomeIntervals data = new GenomeIntervals(List.of(intervals1, intervals2, intervals3));
        List<ChrBaseRegion> expected = new ArrayList<>();
        expected.addAll(regions1);
        expected.addAll(regions2);
        assertEquals(expected, data.regionsList());
    }
}
