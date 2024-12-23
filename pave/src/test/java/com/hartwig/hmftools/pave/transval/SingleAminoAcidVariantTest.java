package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import org.junit.Assert;
import org.junit.Test;

public class SingleAminoAcidVariantTest extends TransvalTestBase
{

    @Test
    public void referenceAminoAcidTest()
    {
        SingleAminoAcidVariant v1 = this.variant("BRAF:p.1E");
        Assert.assertEquals("M", v1.referenceAminoAcid());

        SingleAminoAcidVariant v2 = this.variant("BRAF:p.2E");
        Assert.assertEquals("A", v2.referenceAminoAcid());

        SingleAminoAcidVariant v600 = this.variant("BRAF:p.600E");
        Assert.assertEquals("V", v600.referenceAminoAcid());

        SingleAminoAcidVariant v765 = this.variant("BRAF:p.765E");
        Assert.assertEquals("V", v765.referenceAminoAcid());

        SingleAminoAcidVariant v766 = this.variant("BRAF:p.766E");
        Assert.assertEquals("H", v766.referenceAminoAcid());
    }

    @Test
    public void variantAminoAcidTest()
    {
        SingleAminoAcidVariant v1 = this.variant("BRAF:p.1E");
        Assert.assertEquals("E", v1.variantAminoAcid());
    }

    @Test
    public void possibleVariantCodonsTest()
    {
        SingleAminoAcidVariant v1 = this.variant("BRAF:p.1F");
        Assert.assertEquals(Set.of("TTC", "TTT"), v1.possibleVariantCodons());

        SingleAminoAcidVariant v7 = this.variant("BRAF:p.7G");
        Assert.assertEquals(Set.of("GGA", "GGC", "GGG", "GGT"), v7.possibleVariantCodons());
    }
}
