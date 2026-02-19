package com.hartwig.hmftools.amber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;

import org.junit.Test;

public class AmberSiteTest
{
    @Test
    public void rawPositionTest()
    {
        AmberSite site = new AmberSite("chr1", 100, "A", "G", false);
        GenomePositionImpl rawPosition = site.rawPosition();
        assertEquals("chr1", rawPosition.chromosome());
        assertEquals(100, rawPosition.position());
    }

    @Test
    public void rawPositionsConstantTest()
    {
        AmberSite site1 = new AmberSite("1", 1000, "A", "G", false);
        AmberSite site2 = new AmberSite("1", 1000, "A", "C", false);
        assertEquals(site1.rawPosition(), site1.rawPosition());
        assertEquals(site1.rawPosition(), site2.rawPosition());
    }
}
