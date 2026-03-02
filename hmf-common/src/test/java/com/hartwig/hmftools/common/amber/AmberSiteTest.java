package com.hartwig.hmftools.common.amber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.Test;

public class AmberSiteTest
{
    @Test
    public void defaultFrequencyTest()
    {
        AmberSite site = new AmberSite("chr1", 100, "A", "T", true);
        assertEquals(0.5, site.VariantAlleleFrequency, 0.00001);
    }

    @Test
    public void customFrequencyTest()
    {
        AmberSite site = new AmberSite("chr1", 100, "A", "T", true, 0.1);
        assertEquals(0.1, site.VariantAlleleFrequency, 0.00001);
    }

    @Test
    public void equalsTest()
    {
        AmberSite site1 = new AmberSite("chr1", 100, "A", "T", true);
        AmberSite site2 = new AmberSite("chr1", 100, "A", "T", true);
        AmberSite site3 = new AmberSite("chr1", 100, "A", "T", true, 0.49);
        AmberSite site4 = new AmberSite("chr1", 100, "A", "T", false, 0.49);
        AmberSite site5 = new AmberSite("chr1", 100, "A", "G", true, 0.49);
        AmberSite site6 = new AmberSite("chr1", 100, "C", "T", true, 0.49);
        AmberSite site7 = new AmberSite("chr1", 101, "A", "T", true, 0.49);
        AmberSite site8 = new AmberSite("chr2", 100, "A", "T", true, 0.49);
        assertEquals(site1, site2);
        assertEquals(site1, site3);
        assertEquals(site1, site4);
        assertNotEquals(site1, site5);
        assertNotEquals(site1, site6);
        assertNotEquals(site1, site7);
        assertNotEquals(site1, site8);
    }

    @Test
    public void hashCodeTest()
    {
        AmberSite site1 = new AmberSite("chr1", 100, "A", "T", true, 0.6);
        AmberSite site2 = new AmberSite("chr1", 100, "A", "T", true, 0.6);
        assertEquals(site1.hashCode(), site2.hashCode());
    }

    @Test
    public void toStringTest()
    {
        AmberSite site = new AmberSite("chr1", 100, "A", "T", true, 0.61274872);
        assertEquals("chr1:100 A>T snpcheck 0.613", site.toString());
    }
}
