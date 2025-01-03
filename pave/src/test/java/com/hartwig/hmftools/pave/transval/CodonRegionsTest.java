package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Assert;
import org.junit.Test;

public class CodonRegionsTest
{
    private final RefGenomeInterface genome1 = new FixedStringGenome("AAAACCCGGTGGGGGGGGGGATGGGGGGGGGGGGG");
    private final String chr1 = "chr1";

    @Test
    public void retrieveCodonExonsTest()
    {
        CodonRegions cd1 = new CodonRegions(5, ed(5, 10), ed(20, 28, 2));
        Assert.assertEquals("CCG", cd1.retrieveCodon(genome1));

        cd1 = new CodonRegions(5, ed(5, 10), null);
        Assert.assertEquals("CCG", cd1.retrieveCodon(genome1));

        CodonRegions cd2 = new CodonRegions(6, ed(4, 11), null);
        Assert.assertEquals("CGG", cd2.retrieveCodon(genome1));

        CodonRegions cd3 = new CodonRegions(8, ed(5, 10), null);
        Assert.assertEquals("GTG", cd3.retrieveCodon(genome1));

        CodonRegions cd4 = new CodonRegions(9, ed(5, 10), ed(20, 28, 2));
        Assert.assertEquals("TGA", cd4.retrieveCodon(genome1));

        CodonRegions cd5 = new CodonRegions(10, ed(5, 10), ed(20, 28, 2));
        Assert.assertEquals("GAT", cd5.retrieveCodon(genome1));
    }

    private ChrBaseRegion ed(int start, int end, int rank)
    {
        return new ChrBaseRegion(chr1, start, end);
    }

    private ChrBaseRegion ed(int start, int end)
    {
        return ed(start, end, 1);
    }
}
