package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;

import org.junit.Before;
import org.junit.Test;

public class ScaleFusionTest
{
    private List<GenomeRegion> exons = Lists.newArrayList();
    private ScaleFusion victim;

    @Before
    public void setup()
    {
        exons.add(create(100, 2000));
        exons.add(create(3001, 4600));
        exons.add(create(5601, 11000));
        victim = new ScaleFusion(exons);
    }

    @Test
    public void testIntrons()
    {
        final List<GenomeRegion> introns = ScaleFusion.introns(exons);
        assertEquals(2, introns.size());
        assertEquals(create(2001, 3000), introns.get(0));
        assertEquals(create(4601, 5600), introns.get(1));
    }

    @Test
    public void testScaleExons()
    {
        GenomeRegion scaled;

        scaled = victim.scaled(exons.get(0));
        assertEquals(exons.get(0), scaled);

        scaled = victim.scaled(exons.get(1));
        assertEquals(create(2101, 3700), scaled);

        scaled = victim.scaled(exons.get(2));
        assertEquals(create(3801, 9200), scaled);
    }

    @Test
    public void testScaleProteinDomainInOneExon()
    {
        GenomeRegion proteinDomain = create(3501, 4100);
        GenomeRegion scaled = victim.scaled(proteinDomain);
        assertEquals(create(2601, 3200), scaled);
    }

    @Test
    public void testScaleProteinDomainOverTwoExon()
    {
        GenomeRegion proteinDomain = create(1000, 4100);
        GenomeRegion scaled = victim.scaled(proteinDomain);
        assertEquals(create(1000, 3200), scaled);
    }

    @Test
    public void testScaleProteinDomainStartInIntron()
    {
        GenomeRegion proteinDomain = create(2750, 4000);
        GenomeRegion scaled = victim.scaled(proteinDomain);
        assertEquals(create(2075, 3100), scaled);
    }

    @Test
    public void testScaleProteinDomainEndInIntron()
    {
        GenomeRegion proteinDomain = create(3001, 4850);
        GenomeRegion scaled = victim.scaled(proteinDomain);
        assertEquals(create(2101, 3725), scaled);
    }

    @Test
    public void testScaleProteinDomainStartAndEndInIntron()
    {
        GenomeRegion proteinDomain = create(2750, 4850);
        GenomeRegion scaled = victim.scaled(proteinDomain);
        assertEquals(create(2075, 3725), scaled);
    }

    private static GenomeRegion create(long start, long end)
    {
        return GenomeRegions.create("1", start, end);
    }

}
