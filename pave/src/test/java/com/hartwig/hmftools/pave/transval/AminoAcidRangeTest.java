package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Assert;
import org.junit.Test;

public class AminoAcidRangeTest extends TransvalTestBase
{
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
        Assert.assertTrue(aar(1, "M", 2, "A").applies(brafCanonical));
        Assert.assertTrue(aar(1, "M", 6, "G").applies(brafCanonical));
        Assert.assertTrue(aar(4, "L", 6, "G").applies(brafCanonical));
        Assert.assertTrue(aar(4, "L", 6, null).applies(brafCanonical));
        Assert.assertTrue(aar(4, null, 6, "G").applies(brafCanonical));
        Assert.assertTrue(aar(4, null, 6, null).applies(brafCanonical));

        Assert.assertTrue(aar(4, "L", brafLength - 1, "H").applies(brafCanonical));
        Assert.assertTrue(aar(4, "L", brafLength, "X").applies(brafCanonical));
        Assert.assertFalse(aar(4, "L", brafLength + 1, "H").applies(brafCanonical));
    }

    @Test
    public void startPositionTest()
    {
        Assert.assertEquals(1, aar(1, "M", 2, "A").startPosition());
        Assert.assertEquals(4, aar(4, null, 6, null).startPosition());
    }

    @Test
    public void aminoAcidAtStartTest()
    {
        Assert.assertEquals(aa("M"), aar(1, "M", 2, "A").aminoAcidAtStart());
    }

    @Test
    public void lengthTest()
    {
        Assert.assertEquals(2, aar(1, "M", 2, "A").length());
        Assert.assertEquals(5, aar(4, null, 8, "A").length());
    }

    @Test
    public void equalsTest()
    {
        Assert.assertEquals(aar(1, "M", 2, "A"), aar(1, "M", 2, "A"));
        Assert.assertEquals(aar(1, "M", 6, "G"), aar(1, "M", 6, "G"));
        Assert.assertEquals(aar(1, null, 6, "G"), aar(1, null, 6, "G"));
        Assert.assertEquals(aar(1, "M", 6, null), aar(1, "M", 6, null));
        Assert.assertEquals(aar(1, null, 6, null), aar(1, null, 6, null));

        Assert.assertNotEquals(aar(1, "M", 3, "A"), aar(1, "M", 2, "A"));
        Assert.assertNotEquals(aar(1, "M", 2, "F"), aar(1, "M", 2, "A"));
        Assert.assertNotEquals(aar(1, "M", 3, "F"), aar(1, "M", 2, "A"));
        Assert.assertNotEquals(aar(1, "M", 3, "A"), aar(1, "M", 2, null));
        Assert.assertNotEquals(aar(1, "M", 3, null), aar(1, "M", 2, "A"));
    }

    @Test
    public void hashCodeTest()
    {
        Assert.assertEquals(aar(1, "M", 6, "G").hashCode(), aar(1, "M", 6, "G").hashCode());
        Assert.assertEquals(aar(1, null, 6, null).hashCode(), aar(1, null, 6, null).hashCode());
    }

    @Test
    public void toStringTest()
    {
        Assert.assertEquals("[1,M]_[6,G]", aar(1, "M", 6, "G").toString());
        Assert.assertEquals("[1,?]_[6,?]", aar(1, null, 6, null).toString());
    }

    private AminoAcidRange aar(int start, String left, int end, String right)
    {
        return new AminoAcidRange(aas(start, left), aas(end, right));
    }
}
