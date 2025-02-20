package com.hartwig.hmftools.pave.transval;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class SingleAminoAcidVariantTest extends VariantTest
{

    @Test
    public void variantSequenceTest()
    {
        SingleAminoAcidVariant variant = new SingleAminoAcidVariant(gene, transcript, taa, aas(5, "W"));
        Assert.assertEquals("MAAQWAPAAS", variant.variantSequence().sequence());
    }

    @Test
    public void referenceValueTest()
    {
        // Forward strand examples:
        Assert.assertEquals("M", this.saav("ADCK2:p.1E").referenceAminoAcids());
        Assert.assertEquals("V", this.saav("ADCK2:p.2E").referenceAminoAcids());
        Assert.assertEquals("R", this.saav("ADCK2:p.71E").referenceAminoAcids());
        Assert.assertEquals("K", this.saav("ADCK2:p.421E").referenceAminoAcids());
        Assert.assertEquals("P", this.saav("ADCK2:p.625E").referenceAminoAcids());
        Assert.assertEquals("P", this.saav("ADCK2:p.626E").referenceAminoAcids());

        // Reverse strand examples:
        Assert.assertEquals("M", this.saav("BRAF:p.1E").referenceAminoAcids());
        Assert.assertEquals("A", this.saav("BRAF:p.2E").referenceAminoAcids());
        Assert.assertEquals("V", this.saav("BRAF:p.600E").referenceAminoAcids());
        Assert.assertEquals("V", this.saav("BRAF:p.765E").referenceAminoAcids());
        Assert.assertEquals("H", this.saav("BRAF:p.766E").referenceAminoAcids());
    }

    @Test
    public void codingRegionLengthsTest()
    {
        // TODO cases where there are 1 or 2 exons
        SingleAminoAcidVariant variant = this.saav("ADCK2:p.3K");
        Assert.assertEquals(627, variant.mAminoAcidSequence.AminoAcids.length()); // sanity
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(8, returned.size());
        Assert.assertEquals(3 * 627, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(933, returned.get(0).intValue());
        Assert.assertEquals(147, returned.get(1).intValue());
        Assert.assertEquals(54, returned.get(6).intValue());
        Assert.assertEquals(141, returned.get(7).intValue());
    }

    @Test
    public void codingRegionExons()
    {
        SingleAminoAcidVariant variant = this.saav("ZYX:p.1E");
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(9, returned.size());
        Assert.assertEquals(573, variant.mAminoAcidSequence.AminoAcids.length()); // sanity
        Assert.assertEquals(3 * 573, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(208, returned.get(0).intValue()); // coding start > end of this exon
        Assert.assertEquals(200, returned.get(1).intValue());
        Assert.assertEquals(121, returned.get(7).intValue());
        Assert.assertEquals(105, returned.get(8).intValue());
    }

    @Test
    public void codingRegionLengthsReverseStrandTest()
    {
        SingleAminoAcidVariant variant = this.saav("BRAF:p.600E");
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(18, returned.size());
        Assert.assertEquals(3 * 767, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(138, returned.get(0).intValue());
        Assert.assertEquals(102, returned.get(1).intValue());
        Assert.assertEquals(135, returned.get(16).intValue());
        Assert.assertEquals(174, returned.get(17).intValue());
    }

    @Test
    public void altValueTest()
    {
        Assert.assertEquals("E", this.saav("BRAF:p.1E").altValue());
    }
    @Test
    public void positionOfFirstAlteredCodonTest()
    {
        Assert.assertEquals(1, saav("ZYX:p.1E").positionOfFirstAlteredCodon());
        Assert.assertEquals(2230, saav("MTOR:p.L2230V").positionOfFirstAlteredCodon());
    }

    @Test
    public void positionOfLastAlteredCodonTest()
    {
        Assert.assertEquals(1, saav("ZYX:p.1E").positionOfLastAlteredCodon());
        Assert.assertEquals(2230, saav("MTOR:p.L2230V").positionOfLastAlteredCodon());
    }
}
