package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class AminoAcidTest
{
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
