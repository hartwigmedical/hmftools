package com.hartwig.hmftools.pave.reverse;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

public class SingleAminoAcidVariantTest extends VariantTest
{

    @Test
    public void variantSequenceTest()
    {
        SingleAminoAcidVariant variant = new SingleAminoAcidVariant(gene, transcript, taa, aas(5, "W"));
        assertEquals("MAAQWAPAAS", variant.variantSequence().sequence());
    }

    @Test
    public void codingRegionLengthsTest()
    {
        SingleAminoAcidVariant variant = this.saav("ADCK2:p.V2K");
        assertEquals(627, variant.mAminoAcidSequence.AminoAcids.length()); // sanity
        List<Integer> returned = variant.codingRegionLengths();
        assertEquals(8, returned.size());
        assertEquals(3 * 627, returned.stream().mapToInt(Integer::intValue).sum());
        assertEquals(933, returned.get(0).intValue());
        assertEquals(147, returned.get(1).intValue());
        assertEquals(54, returned.get(6).intValue());
        assertEquals(141, returned.get(7).intValue());
    }

    @Test
    public void codingRegionExons()
    {
        SingleAminoAcidVariant variant = this.saav("ZYX:p.A2E");
        List<Integer> returned = variant.codingRegionLengths();
        assertEquals(9, returned.size());
        assertEquals(573, variant.mAminoAcidSequence.AminoAcids.length()); // sanity
        assertEquals(3 * 573, returned.stream().mapToInt(Integer::intValue).sum());
        assertEquals(208, returned.get(0).intValue()); // coding start > end of this exon
        assertEquals(200, returned.get(1).intValue());
        assertEquals(121, returned.get(7).intValue());
        assertEquals(105, returned.get(8).intValue());
    }

    @Test
    public void codingRegionLengthsReverseStrandTest()
    {
        SingleAminoAcidVariant variant = this.saav("BRAF:p.V600E");
        List<Integer> returned = variant.codingRegionLengths();
        assertEquals(18, returned.size());
        assertEquals(3 * 767, returned.stream().mapToInt(Integer::intValue).sum());
        assertEquals(138, returned.get(0).intValue());
        assertEquals(102, returned.get(1).intValue());
        assertEquals(135, returned.get(16).intValue());
        assertEquals(174, returned.get(17).intValue());
    }

    @Test
    public void altValueTest()
    {
        assertEquals("E", this.saav("VHL:p.P2E").altValue());
    }

    @Test
    public void positionOfFirstAlteredCodonTest()
    {
        assertEquals(2, saav("VHL:p.P2E").positionOfFirstAlteredCodon());
        assertEquals(2230, saav("MTOR:p.L2230V").positionOfFirstAlteredCodon());
    }
}
