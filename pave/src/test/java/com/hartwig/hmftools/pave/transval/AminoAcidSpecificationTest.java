package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Assert;
import org.junit.Test;

public class AminoAcidSpecificationTest extends TransvalTestBase
{
    @Test
    public void equalsTest()
    {
        Assert.assertEquals(aas(1, "A"), aas(1, "A"));
        Assert.assertEquals(aas(1, null), aas(1, null));
        Assert.assertEquals(aas(1, "Q"), aas(1, "Q"));
        Assert.assertNotEquals(aas(10, "Q"), aas(1, "Q"));
        Assert.assertNotEquals(aas(10, "Q"), aas(1, null));
        Assert.assertNotEquals(aas(1, "Q"), aas(1, "A"));
        Assert.assertNotEquals(aas(1, "Q"), aas(10, "A"));
    }

    @Test
    public void hashCodeTest()
    {
        Assert.assertEquals(aas(1, "Q").hashCode(), aas(1, "Q").hashCode());
        Assert.assertEquals(aas(100, null).hashCode(), aas(100, null).hashCode());
    }

    @Test
    public void toStringTest()
    {
        Assert.assertEquals("[1,A]", aas(1,"A").toString());
        Assert.assertEquals("[1,?]", aas(1,null).toString());
    }

    @Test
    public void constructFromString()
    {
        Assert.assertEquals(new AminoAcidSpecification(1, (String) null), aas(1, null));

    }

    @Test
    public void filterTranscriptTest()
    {
        /*
        BRAF, ENST00000646891 (canonical):
        1   2   3   4   5   6   ...
        M   A   A   L   S   G   ... P V H X
         */
        final TranscriptData brafCanonicalTranscript = transcript("ENSG00000157764", "ENST00000646891");
        final TranscriptAminoAcids brafCanonical = transval.mTranscriptAminoAcidsMap.get(brafCanonicalTranscript.TransName);
        final int brafLength = brafCanonical.AminoAcids.length();
        Assert.assertTrue(aas(1, "M").applies(brafCanonical));
        Assert.assertTrue(aas(2, "A").applies(brafCanonical));
        Assert.assertTrue(aas(3, "A").applies(brafCanonical));
        Assert.assertFalse(aas(3, "L").applies(brafCanonical));
        Assert.assertTrue(aas(brafLength - 1, "H").applies(brafCanonical));
        Assert.assertTrue(aas(brafLength, "X").applies(brafCanonical));
        Assert.assertFalse(aas(brafLength + 1, "H").applies(brafCanonical));

        Assert.assertTrue(aas(4, null).applies(brafCanonical));
        Assert.assertTrue(aas(brafLength, null).applies(brafCanonical));
        Assert.assertFalse(aas(brafLength + 1, null).applies(brafCanonical));
    }

    private AminoAcidSpecification aas(int position, String symbol)
    {
        final AminoAcid aa = symbol == null ? null : aa(symbol);
        return new AminoAcidSpecification(position, aa);
    }
}
