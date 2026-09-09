package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.AMP;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.DEL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.junit.Test;

public class RegionAmpDelTest extends ToolsTestBase
{
    @Test
    public void constructorTest()
    {
        RegionGeneEvents rge = rge("G", _1, 1000, 2000, AMPLIFICATION);
        RegionAmpDel ampDel = new RegionAmpDel(rge);
        assertEquals(1000, ampDel.start());
        assertEquals(2000, ampDel.end());
    }

    @Test
    public void equalityTest()
    {
        assertNotEquals(rad(_1, 100, 200, AMPLIFICATION), rad(_1, 100, 200, HET_DELETION));
        assertNotEquals(rad(_1, 100, 200, AMPLIFICATION), rad(_1, 100, 200, HOM_DELETION));
        assertNotEquals(rad(_1, 100, 200, AMPLIFICATION), rad(_1, 100, 201, AMPLIFICATION));
        assertNotEquals(rad(_1, 100, 200, AMPLIFICATION), rad(_1, 99, 200, AMPLIFICATION));
        assertNotEquals(rad(_1, 100, 200, AMPLIFICATION), rad(_2, 100, 200, AMPLIFICATION));
        assertNotEquals(rad(_1, 100, 200, AMPLIFICATION), "whatever");
        assertNotEquals(rad(_1, 100, 200, AMPLIFICATION), null);

        assertEquals(rad(_1, 100, 200, AMPLIFICATION), rad(_1, 100, 200, AMPLIFICATION));
        assertEquals(rad(_1, 100, 200, AMPLIFICATION).hashCode(), rad(_1, 100, 200, AMPLIFICATION).hashCode());
    }

    @Test
    public void compareToTest()
    {
        assertEquals(AMP.compareTo(DEL), rad(_1, 100, 200, AMPLIFICATION).compareTo(rad(_1, 100, 200, HET_DELETION)));
        assertEquals(_1.compareTo(_2), rad(_1, 100, 200, AMPLIFICATION).compareTo(rad(_2, 100, 200, AMPLIFICATION)));
        assertEquals(-1, rad(_1, 100, 200, AMPLIFICATION).compareTo(rad(_1, 100, 201, AMPLIFICATION)));
        assertEquals(-1, rad(_1, 100, 200, AMPLIFICATION).compareTo(rad(_1, 101, 200, AMPLIFICATION)));
        assertEquals(0, rad(_1, 100, 200, AMPLIFICATION).compareTo(rad(_1, 100, 200, AMPLIFICATION)));
    }

    @Test
    public void typeTest()
    {
        assertEquals("AMP", rad(_1, 100, 200, AMPLIFICATION).type());
        assertEquals("DEL", rad(_1, 100, 200, HET_DELETION).type());
    }

    private RegionAmpDel rad(HumanChromosome chromosome, int start, int end, GermlineStatus germlineStatus)
    {
        return new RegionAmpDel(rge("A", chromosome, start, end, germlineStatus));
    }
}
