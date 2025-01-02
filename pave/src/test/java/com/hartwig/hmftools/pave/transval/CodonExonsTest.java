package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.junit.Assert;
import org.junit.Test;

public class CodonExonsTest
{
    private final RefGenomeInterface genome1 = new FixedStringGenome("AAAACCCGGTGGGGGGGGGGATGGGGGGGGGGGGG");
    private final String chr1 = "chr1";

    @Test
    public void retrieveCodonExonsTest()
    {
        CodonExons cd1 = new CodonExons(5, ed(5, 10), ed(20, 28, 2));
        Assert.assertEquals("CCG", cd1.retrieveCodon(chr1, genome1));

        cd1 = new CodonExons(5, ed(5, 10), null);
        Assert.assertEquals("CCG", cd1.retrieveCodon(chr1, genome1));

        CodonExons cd2 = new CodonExons(6, ed(4, 11), null);
        Assert.assertEquals("CGG", cd2.retrieveCodon(chr1, genome1));

        CodonExons cd3 = new CodonExons(8, ed(5, 10), null);
        Assert.assertEquals("GTG", cd3.retrieveCodon(chr1, genome1));

        CodonExons cd4 = new CodonExons(9, ed(5, 10), ed(20, 28, 2));
        Assert.assertEquals("TGA", cd4.retrieveCodon(chr1, genome1));

        CodonExons cd5 = new CodonExons(10, ed(5, 10), ed(20, 28, 2));
        Assert.assertEquals("GAT", cd5.retrieveCodon(chr1, genome1));
    }

    private ExonData ed(int start, int end, int rank)
    {
        return new ExonData(999, start, end, rank, start, end);
    }

    private ExonData ed(int start, int end)
    {
        return ed(start, end, 1);
    }
}
