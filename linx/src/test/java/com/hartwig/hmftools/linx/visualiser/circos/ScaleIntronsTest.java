package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.junit.Before;
import org.junit.Test;

public class ScaleIntronsTest
{
    private List<GenomeRegion> exons = Lists.newArrayList();
    private ScaleIntrons victim;

    private int[] original = {100, 2000, 3001, 4600, 5601, 11000};
    private int[] scaled = {100, 2000, 2011, 3610, 3621, 9020};

    @Before
    public void setup()
    {
        exons.add(create(original[0], original[1]));
        exons.add(create(original[2], original[3]));
        exons.add(create(original[4], original[5]));
        victim = new ScaleIntrons(ScaleIntrons.introns(exons));
    }

    @Test
    public void testIntrons()
    {
        final List<GenomeRegion> introns = ScaleIntrons.introns(exons);
        assertEquals(2, introns.size());
        assertEquals(create(2001, 3000), introns.get(0));
        assertEquals(create(4601, 5600), introns.get(1));
    }

    @Test
    public void testScaleExons()
    {
        GenomeRegion region;

        region = victim.scaleIntronsFromRegion(exons.get(0));
        assertEquals(create(scaled[0], scaled[1]), region);

        region = victim.scaleIntronsFromRegion(exons.get(1));
        assertEquals(create(scaled[2], scaled[3]), region);

        region = victim.scaleIntronsFromRegion(exons.get(2));
        assertEquals(create(scaled[4], scaled[5]), region);
    }

    @Test
    public void testScaleProteinDomainInOneExon()
    {
        GenomeRegion proteinDomain = create(original[2] + 10, original[3] - 10);
        GenomeRegion region = victim.scaleIntronsFromRegion(proteinDomain);
        assertEquals(create(scaled[2] + 10, scaled[3] - 10), region);
    }

    @Test
    public void testScaleProteinDomainOverTwoExon()
    {
        GenomeRegion proteinDomain = create(original[0] + 10, original[3] - 10);
        GenomeRegion region = victim.scaleIntronsFromRegion(proteinDomain);
        assertEquals(create(scaled[0] + 10, scaled[3] - 10), region);
    }

    @Test
    public void testScaleProteinDomainStartInIntron()
    {
        GenomeRegion proteinDomain = create(original[2] - 250, original[3] - 10);
        GenomeRegion region = victim.scaleIntronsFromRegion(proteinDomain);
        assertEquals(create(scaled[2] - 3, scaled[3] - 10), region);
    }

    @Test
    public void testScaleProteinDomainEndInIntron()
    {
        GenomeRegion proteinDomain = create(original[2] + 10, original[3] + 250);
        GenomeRegion region = victim.scaleIntronsFromRegion(proteinDomain);
        assertEquals(create(scaled[2] + 10, scaled[3] + 3), region);
    }

    @Test
    public void testScaleProteinDomainStartAndEndInIntron()
    {
        GenomeRegion proteinDomain = create(original[2] - 250, original[3] + 250);
        GenomeRegion region = victim.scaleIntronsFromRegion(proteinDomain);
        assertEquals(create(scaled[2] - 3, scaled[3]  + 3), region);
    }

    private static GenomeRegion create(long start, long end)
    {
        return GenomeRegions.create("1", start, end);
    }

}
