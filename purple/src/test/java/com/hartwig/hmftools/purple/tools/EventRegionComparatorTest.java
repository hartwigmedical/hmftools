package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class EventRegionComparatorTest extends ToolsTestBase
{
    EventRegionComparator comparator = new EventRegionComparator();

    @Test
    public void compareToTest()
    {
        assertEquals(_1.compareTo(_2), comparator.compare(cbr(_1, 100, 200), cbr(_2, 100, 200)));
        assertEquals(_3.compareTo(_2), comparator.compare(cbr(_3, 100, 200), cbr(_2, 100, 200)));
        assertEquals(Integer.compare(100, 101), comparator.compare(cbr(_3, 100, 200), cbr(_3, 101, 200)));
        assertEquals(Integer.compare(101, 100), comparator.compare(cbr(_3, 101, 200), cbr(_3, 100, 200)));
        assertEquals(Integer.compare(200, 201), comparator.compare(cbr(_3, 100, 200), cbr(_3, 100, 201)));
        assertEquals(Integer.compare(201, 200), comparator.compare(cbr(_3, 100, 201), cbr(_3, 100, 200)));

        assertEquals(0, comparator.compare(cbr(_3, 100, 200), cbr(_3, 100, 200)));
    }
}
