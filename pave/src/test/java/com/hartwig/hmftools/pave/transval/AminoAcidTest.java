package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class AminoAcidTest extends TransvalTestBase
{
    private final AminoAcid A = aa("A"); //GCT, GCA, GCC, GCG
    private final AminoAcid R = aa("R"); // AGG, AGA, CGA, CGC, CGG, CGT
    private final AminoAcid W = aa("W"); // TGG

    @Test
    public void matchingTruncatedCodons()
    {
        Assert.assertEquals(css("GCT, GCA, GCC, GCG"), A.matchingTruncatedCodons("", ""));
        Assert.assertEquals(css("CT, CA, CC, CG"), A.matchingTruncatedCodons("G", ""));
        Assert.assertEquals(css("T, A, C, G"), A.matchingTruncatedCodons("GC", ""));
        Assert.assertTrue(A.matchingTruncatedCodons("GC", "G").isEmpty());
        Assert.assertEquals(css("C"), A.matchingTruncatedCodons("G", "G"));
        Assert.assertEquals(css("GC"), A.matchingTruncatedCodons("", "G"));
        Assert.assertEquals(css("G"), A.matchingTruncatedCodons("", "CG"));
        Assert.assertTrue(A.matchingTruncatedCodons("T", "").isEmpty());

        Assert.assertEquals(css("AGG, AGA, CGA, CGC, CGG, CGT"), R.matchingTruncatedCodons("", ""));
        Assert.assertEquals(css("GG, GA"), R.matchingTruncatedCodons("A", ""));
        Assert.assertEquals(css("G, A"), R.matchingTruncatedCodons("AG", ""));

        Assert.assertEquals(css("TGG"), W.matchingTruncatedCodons("", ""));
        Assert.assertEquals(css("TG"), W.matchingTruncatedCodons("", "G"));
    }

    @Test
    public void construction()
    {
        checkError("a");
        checkError("junk");
        checkError("Alanine");
        checkError("AA");
        checkError("Z"); // means 'Glutamine or Glutamic Acid'
        checkError("Glx"); // means 'Glutamine or Glutamic Acid'
        checkError("B"); //// B and Asx mean 'Asparagine or Aspartic Acid'
        checkError("Asx"); //// B and Asx mean 'Asparagine or Aspartic Acid'
    }

    @Test
    public void isConvertedToSingleLetter()
    {
        Assert.assertEquals("A", new AminoAcid("A").symbol);
        Assert.assertEquals("A", new AminoAcid("Ala").symbol);
    }

    @Test
    public void equalsTest()
    {
        Assert.assertEquals(new AminoAcid("A"), new AminoAcid("A"));
        Assert.assertEquals(new AminoAcid("A"), new AminoAcid("Ala"));
        Assert.assertEquals(new AminoAcid("Thr"), new AminoAcid("T"));
        Assert.assertNotEquals(new AminoAcid("A"), new AminoAcid("T"));
    }

    @Test
    public void hashCodeTest()
    {
        Assert.assertEquals(new AminoAcid("A").hashCode(), new AminoAcid("A").hashCode());
        Assert.assertEquals(new AminoAcid("W").hashCode(), new AminoAcid("Trp").hashCode());
    }

    @Test
    public void toStringTest()
    {
        Assert.assertEquals("A", new AminoAcid("A").toString());
    }

    private void checkError(String s)
    {
        try
        {
            new AminoAcid(s);
            Assert.fail("Error should have been thrown for '" + s + "'");
        }
        catch(final Exception e)
        {
            // ignore
        }
    }
}
