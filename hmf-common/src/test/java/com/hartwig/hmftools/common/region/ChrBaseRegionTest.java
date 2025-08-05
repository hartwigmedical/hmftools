package com.hartwig.hmftools.common.region;

import static java.util.List.of;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Assert;
import org.junit.Test;

public class ChrBaseRegionTest
{
    private final ChrBaseRegion theRegion = new ChrBaseRegion("12", 1000, 2000);

    @Test
    public void splitBySingleOverlapTest()
    {
        ChrBaseRegion gr0 = cbr(900, 2100, theRegion.chromosome());

        List<Pair<ChrBaseRegion, ChrBaseRegion>> found = theRegion.splitByOverlappingRegions(of(gr0));
        assertEquals(1, found.size());
        assertEquals(theRegion, found.get(0).getKey());
        assertEquals(gr0, found.get(0).getValue());
    }

    @Test
    public void splitByOverlapsTest()
    {
        ChrBaseRegion gr1 = cbr(951, 1050, theRegion.chromosome());
        ChrBaseRegion gr2 = cbr(1051, 1250, theRegion.chromosome());
        ChrBaseRegion gr3 = cbr(1251, 1550, theRegion.chromosome());
        ChrBaseRegion gr4 = cbr(1551, 2500, theRegion.chromosome());

        List<Pair<ChrBaseRegion, ChrBaseRegion>> found = theRegion.splitByOverlappingRegions(of(gr1, gr2, gr3, gr4));
        assertEquals(4, found.size());
        assertEquals(Pair.of(cbr(theRegion.start(), gr1.end(), theRegion.Chromosome), gr1), found.get(0));
        assertEquals(Pair.of(gr2, gr2), found.get(1));
        assertEquals(Pair.of(gr3, gr3), found.get(2));
        assertEquals(Pair.of(cbr(gr4.start(), theRegion.end(), theRegion.Chromosome), gr4), found.get(3));
    }

    @Test
    public void findIntersectingRegionsWithEmpty()
    {
        Assert.assertTrue(theRegion.findOverlaps(of()).isEmpty());
    }

    @Test
    public void findIntersectingRegionsDifferentChromosome()
    {
        Assert.assertTrue(theRegion.findOverlaps(of(cbr(900, 1100, "chrZ"))).isEmpty());
    }

    @Test
    public void findIntersectingRegionsTotalCoverageByASingleRegion()
    {
        ChrBaseRegion gr0 = cbr(900, 950, theRegion.chromosome());
        ChrBaseRegion gr1 = cbr(951, 2500, theRegion.chromosome());
        ChrBaseRegion gr2 = cbr(2501, 5500, theRegion.chromosome());

        List<ChrBaseRegion> found = theRegion.findOverlaps(of(gr0, gr1, gr2));
        assertEquals(1, found.size());
        assertEquals(gr1, found.get(0));
    }

    @Test
    public void findIntersectingRegionsFirstRegionEndsAtStart()
    {
        ChrBaseRegion gr0 = cbr(900, 1000, theRegion.chromosome());
        ChrBaseRegion gr1 = cbr(1001, 2500, theRegion.chromosome());
        ChrBaseRegion gr2 = cbr(2501, 5500, theRegion.chromosome());

        List<ChrBaseRegion> found = theRegion.findOverlaps(of(gr0, gr1, gr2));
        assertEquals(2, found.size());
        assertEquals(gr0, found.get(0));
        assertEquals(gr1, found.get(1));
    }

    @Test
    public void findIntersectingRegionsFirstRegionEndsAtEnd()
    {
        ChrBaseRegion gr0 = cbr(900, 2000, theRegion.chromosome());
        ChrBaseRegion gr1 = cbr(2001, 2500, theRegion.chromosome());
        ChrBaseRegion gr2 = cbr(2501, 5500, theRegion.chromosome());

        List<ChrBaseRegion> found = theRegion.findOverlaps(of(gr0, gr1, gr2));
        assertEquals(1, found.size());
        assertEquals(gr0, found.get(0));
    }

    @Test
    public void severalRegionsIntersect()
    {
        ChrBaseRegion gr0 = cbr(900, 950, theRegion.chromosome());
        ChrBaseRegion gr1 = cbr(951, 1050, theRegion.chromosome());
        ChrBaseRegion gr2 = cbr(1051, 1250, theRegion.chromosome());
        ChrBaseRegion gr3 = cbr(1251, 1550, theRegion.chromosome());
        ChrBaseRegion gr4 = cbr(1551, 2500, theRegion.chromosome());
        ChrBaseRegion gr5 = cbr(2551, 3500, theRegion.chromosome());

        List<ChrBaseRegion> found = theRegion.findOverlaps(of(gr0, gr1, gr2, gr3, gr4, gr5));
        assertEquals(4, found.size());
        assertEquals(gr1, found.get(0));
        assertEquals(gr2, found.get(1));
        assertEquals(gr3, found.get(2));
        assertEquals(gr4, found.get(3));
    }

    @Test
    public void overlapsTest()
    {
        assertFalse(theRegion.overlaps(cbr(900, 950, theRegion.chromosome())));
        assertFalse(cbr(900, 950, theRegion.chromosome()).overlaps(theRegion));
        assertFalse(theRegion.overlaps(cbr(2050, 2250, theRegion.chromosome())));
        assertFalse(cbr(2050, 2250, theRegion.chromosome()).overlaps(theRegion));

        checkIntersects(theRegion, cbr(950, 1000, theRegion.chromosome()));
        checkIntersects(theRegion, cbr(950, 1050, theRegion.chromosome()));
        checkIntersects(theRegion, cbr(1050, 1950, theRegion.chromosome()));
        checkIntersects(theRegion, cbr(1950, 2050, theRegion.chromosome()));
        checkIntersects(theRegion, cbr(1999, 2050, theRegion.chromosome()));
        checkIntersects(theRegion, cbr(2000, 2050, theRegion.chromosome()));

        assertFalse(theRegion.overlaps(cbr(2001, 2050, theRegion.chromosome())));

        assertFalse(theRegion.overlaps(cbr(1050, 1950, "ChrU")));
    }

    private ChrBaseRegion cbr(int start, int end, String chromosome)
    {
        return new ChrBaseRegion(chromosome, start, end);
    }

    private static void checkIntersects(ChrBaseRegion region1, ChrBaseRegion region2)
    {
        assertTrue(region1.overlaps(region2));
        assertTrue(region2.overlaps(region1));
    }
}
