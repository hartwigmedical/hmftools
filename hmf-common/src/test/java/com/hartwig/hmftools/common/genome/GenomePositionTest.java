package com.hartwig.hmftools.common.genome;

import java.util.List;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.junit.Assert;
import org.junit.Test;

public class GenomePositionTest
{
    private final GenomePosition theSpot = new GP(1000, "chr1");

    @Test
    public void filterEmpty()
    {
        Assert.assertTrue(theSpot.findContainingRegions(List.of()).isEmpty());
    }

    @Test
    public void differentChromosome()
    {
        Assert.assertTrue(theSpot.findContainingRegions(List.of(new GR(900, 1100, "chrZ"))).isEmpty());
    }

    @Test
    public void atStart()
    {
        List<GenomeRegion> found = theSpot.findContainingRegions(List.of(new GR(1000, 1100, theSpot.chromosome())));
        Assert.assertEquals(1, found.size());
    }

    @Test
    public void atEnd()
    {
        List<GenomeRegion> found = theSpot.findContainingRegions(List.of(new GR(900, 1000, theSpot.chromosome())));
        Assert.assertEquals(1, found.size());
    }

    @Test
    public void mixed()
    {
        GR gr0 = new GR(900, 950, theSpot.chromosome());
        GR gr1 = new GR(951, 1500, theSpot.chromosome());
        GR gr2 = new GR(961, 1500, theSpot.chromosome());
        GR gr3 = new GR(971, 1500, theSpot.chromosome());
        GR gr4 = new GR(1501, 2000, theSpot.chromosome());
        GR gr5 = new GR(951, 1500, "X");
        GR gr6 = new GR(951, 1500, "Y");

        List<GenomeRegion> found = theSpot.findContainingRegions(List.of(gr0, gr1, gr2, gr3, gr4, gr5, gr6));
        Assert.assertEquals(3, found.size());
        Assert.assertEquals(gr1, found.get(0));
        Assert.assertEquals(gr2, found.get(1));
        Assert.assertEquals(gr3, found.get(2));
    }
}

record GP(int position, String chromosome) implements GenomePosition
{

}

record GR(int start, int end, String chromosome) implements GenomeRegion
{

}
