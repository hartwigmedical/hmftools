package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.AMP;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.DEL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class RegionGeneEventsTest extends ToolsTestBase
{
    @Test
    public void constructorTest()
    {
        RegionGeneEvents rge = rge("G", _1, 1000, 2000, AMPLIFICATION);
        assertEquals(rge.genes(), List.of("G"));
        assertEquals(_1, rge.chromosome());
        assertEquals(1000, rge.region().start());
        assertEquals(2000, rge.region().end());
        assertEquals(_1, rge.region().humanChromosome());
        assertEquals(AMP, rge.eventType());
    }

    @Test
    public void eventTypeTest()
    {
        RegionGeneEvents rge = rge("G", _1, 1000, 2000, HET_DELETION);
        assertEquals(DEL, rge.eventType());

        rge = rge("G", _1, 1000, 2000, HOM_DELETION);
        assertEquals(DEL, rge.eventType());
    }

    @Test
    public void offerRejectedTest()
    {
        RegionGeneEvents rge1 = rge("G", _1, 1000, 2000, HET_DELETION);
        assertFalse(rge1.offer(gad("G", _2, 1000, 2000, HET_DELETION)));
        assertEquals(1, rge1.genes().size());
        assertFalse(rge1.offer(gad("G", _1, 1100, 2000, HET_DELETION)));
        assertEquals(1, rge1.genes().size());
        assertFalse(rge1.offer(gad("G", _1, 1000, 2001, HET_DELETION)));
        assertEquals(1, rge1.genes().size());
        assertFalse(rge1.offer(gad("G", _1, 1000, 2000, AMPLIFICATION)));
        assertEquals(1, rge1.genes().size());
    }

    @Test
    public void offerAcceptedTest()
    {
        RegionGeneEvents rge1 = rge("G", _1, 1000, 2000, HET_DELETION);
        assertTrue(rge1.offer(gad("H", _1, 1000, 2000, HET_DELETION)));
        assertEquals(2, rge1.genes().size());
        assertTrue(rge1.genes().contains("G"));
        assertTrue(rge1.genes().contains("H"));
    }

    @Test
    public void offerThrowsExceptionIfOrderNotFollowedTest()
    {
        RegionGeneEvents rge1 = rge("G", _2, 1000, 2000, HET_DELETION);
        assertThrows(IllegalArgumentException.class, () -> rge1.offer(gad("H", _1, 1000, 2000, HET_DELETION)));
        assertThrows(IllegalArgumentException.class, () -> rge1.offer(gad("H", _2, 999, 2000, HET_DELETION)));
        assertThrows(IllegalArgumentException.class, () -> rge1.offer(gad("H", _2, 1000, 1999, HET_DELETION)));
        assertThrows(IllegalArgumentException.class, () -> rge1.offer(gad("H", _2, 1000, 2000, HOM_DELETION)));
    }
}
