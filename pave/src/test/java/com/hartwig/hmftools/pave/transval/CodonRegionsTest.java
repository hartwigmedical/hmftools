package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Assert;
import org.junit.Test;

public class CodonRegionsTest
{
    private final RefGenomeInterface genome1 = new FixedStringGenome("AAAACCCGGTGGGGGGGGGGATGGGGGGGGGGGGG");

    @Test
    public void retrieveCodonsPositiveStrandTest()
    {
        CodonRegions cd1 = new CodonRegions(5, cbr(5, 10), cbr(20, 28));
        Assert.assertEquals("CCG", cd1.retrieveCodon(genome1));

        cd1 = new CodonRegions(5, cbr(5, 10), null);
        Assert.assertEquals("CCG", cd1.retrieveCodon(genome1));

        CodonRegions cd2 = new CodonRegions(6, cbr(4, 11), null);
        Assert.assertEquals("CGG", cd2.retrieveCodon(genome1));

        CodonRegions cd3 = new CodonRegions(8, cbr(5, 10), null);
        Assert.assertEquals("GTG", cd3.retrieveCodon(genome1));

        CodonRegions cd4 = new CodonRegions(9, cbr(5, 10), cbr(20, 28));
        Assert.assertEquals("TGA", cd4.retrieveCodon(genome1));

        CodonRegions cd5 = new CodonRegions(10, cbr(5, 10), cbr(20, 28));
        Assert.assertEquals("GAT", cd5.retrieveCodon(genome1));
    }

    @Test
    public void retrieveCodonsNegativeStrandTest()
    {
        // -----CCGGTG---------GATGGGGGG--------
        //      5    10        20     28
        CodonRegions cd1 = new CodonRegions(7, cbr(5, 10), null, false);
        Assert.assertEquals("CGG", cd1.retrieveCodon(genome1));

        CodonRegions cd2 = new CodonRegions(8, cbr(5, 10),null, false);
        Assert.assertEquals("CCG", cd2.retrieveCodon(genome1));

        CodonRegions cd3 = new CodonRegions(9, cbr(5, 10),null, false);
        Assert.assertEquals("ACC", cd3.retrieveCodon(genome1));

        CodonRegions cd4 = new CodonRegions(19, cbr(19, 28), cbr(5, 10), false);
        Assert.assertEquals("CCA", cd4.retrieveCodon(genome1));

        CodonRegions cd5 = new CodonRegions(20, cbr(19, 28), cbr(5, 10), false);
        Assert.assertEquals("TCC", cd5.retrieveCodon(genome1));

        CodonRegions cd6 = new CodonRegions(21, cbr(19, 28), cbr(5, 10), false);
        Assert.assertEquals("ATC", cd6.retrieveCodon(genome1));
    }

    @Test
    public void singleExon()
    {
        Assert.assertTrue( new CodonRegions(5, cbr(5, 10), cbr(20, 28)).codonIsInSingleExon());
        Assert.assertFalse(new CodonRegions(20, cbr(20, 28), cbr(5, 10), false).codonIsInSingleExon());
    }

    @Test
    public void translateCodonPositionTest()
    {
        CodonRegions cd1 = new CodonRegions(7, cbr(5, 10), null, true);
        Assert.assertEquals(7, cd1.translateCodonPosition(0));
        Assert.assertEquals(8, cd1.translateCodonPosition(1));
        Assert.assertEquals(9, cd1.translateCodonPosition(2));

        CodonRegions cd2 = new CodonRegions(8, cbr(5, 10), cbr(15, 20), true);
        Assert.assertEquals(8, cd2.translateCodonPosition(0));
        Assert.assertEquals(9, cd2.translateCodonPosition(1));
        Assert.assertEquals(10, cd2.translateCodonPosition(2));

        CodonRegions cd3 = new CodonRegions(9, cbr(5, 10), cbr(15, 20), true);
        Assert.assertEquals(9, cd3.translateCodonPosition(0));
        Assert.assertEquals(10, cd3.translateCodonPosition(1));
        Assert.assertEquals(15, cd3.translateCodonPosition(2));

        CodonRegions cd4 = new CodonRegions(10, cbr(5, 10), cbr(15, 20), true);
        Assert.assertEquals(10, cd4.translateCodonPosition(0));
        Assert.assertEquals(15, cd4.translateCodonPosition(1));
        Assert.assertEquals(16, cd4.translateCodonPosition(2));
    }

    @Test
    public void translateCodonPositionNegativeStrandTest()
    {
        CodonRegions cd1 = new CodonRegions(7, cbr(5, 10), null, false);
        Assert.assertEquals(7, cd1.translateCodonPosition(0));
        Assert.assertEquals(6, cd1.translateCodonPosition(1));
        Assert.assertEquals(5, cd1.translateCodonPosition(2));

        CodonRegions cd2 = new CodonRegions(8, cbr(5, 10), null, false);
        Assert.assertEquals(8, cd2.translateCodonPosition(0));
        Assert.assertEquals(7, cd2.translateCodonPosition(1));
        Assert.assertEquals(6, cd2.translateCodonPosition(2));

        CodonRegions cd3 = new CodonRegions(15, cbr(15, 20), cbr(5, 10), false);
        Assert.assertEquals(15, cd3.translateCodonPosition(0));
        Assert.assertEquals(10, cd3.translateCodonPosition(1));
        Assert.assertEquals(9, cd3.translateCodonPosition(2));

        CodonRegions cd4 = new CodonRegions(16, cbr(15, 20), cbr(5, 10), false);
        Assert.assertEquals(16, cd4.translateCodonPosition(0));
        Assert.assertEquals(15, cd4.translateCodonPosition(1));
        Assert.assertEquals(10, cd4.translateCodonPosition(2));
    }

    private ChrBaseRegion cbr(int start, int end)
    {
        return new ChrBaseRegion("chr1", start, end);
    }
}
