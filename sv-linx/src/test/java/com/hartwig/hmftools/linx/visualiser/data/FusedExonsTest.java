package com.hartwig.hmftools.linx.visualiser.data;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.junit.Test;

public class FusedExonsTest
{

    @Test
    public void testPositiveStreamBounds() {
        final GenomeRegion gene = GenomeRegions.create("1" , 10000, 15000);
        final GenomeRegion exon =  GenomeRegions.create("1" , 5000, 20000);
        final GenomeRegion expectedEnd =  GenomeRegions.create("1" , 0, 5000);
        assertEquals(expectedEnd, FusedExons.convertRegion(1, gene, exon));
    }

    @Test
    public void testNegativeStreamBounds() {
        final GenomeRegion gene = GenomeRegions.create("1" , 10000, 15000);
        final GenomeRegion exonicBreak =  GenomeRegions.create("1" , 5000, 20000);
        final GenomeRegion expectedEnd =  GenomeRegions.create("1" , 0, 5000);
        assertEquals(expectedEnd, FusedExons.convertRegion(1, gene, exonicBreak));
    }

}
