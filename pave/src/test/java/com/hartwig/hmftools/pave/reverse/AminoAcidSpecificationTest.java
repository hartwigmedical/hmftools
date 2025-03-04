package com.hartwig.hmftools.pave.reverse;

import static com.hartwig.hmftools.pave.reverse.AminoAcidSpecification.parse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Test;

public class AminoAcidSpecificationTest extends ReversePaveTestBase
{
    @Test
    public void parseTest()
    {
        assertEquals(aas(1, "A"), parse("A1"));
        assertEquals(aas(1, "A"), parse("Ala1"));
        assertEquals(aas(1234, "A"), parse("Ala1234"));
        assertEquals(aas(1234, "R"), parse("Arg1234"));
    }

    @Test
    public void equalsTest()
    {
        assertEquals(aas(1, "A"), aas(1, "A"));
        assertEquals(aas(1, null), aas(1, null));
        assertEquals(aas(1, "Q"), aas(1, "Q"));
        assertNotEquals(aas(10, "Q"), aas(1, "Q"));
        assertNotEquals(aas(10, "Q"), aas(1, null));
        assertNotEquals(aas(1, "Q"), aas(1, "A"));
        assertNotEquals(aas(1, "Q"), aas(10, "A"));
    }

    @Test
    public void hashCodeTest()
    {
        assertEquals(aas(1, "Q").hashCode(), aas(1, "Q").hashCode());
        assertEquals(aas(100, null).hashCode(), aas(100, null).hashCode());
    }

    @Test
    public void toStringTest()
    {
        assertEquals("[1,A]", aas(1,"A").toString());
        assertEquals("[1,?]", aas(1,null).toString());
    }

    @Test
    public void constructFromString()
    {
        assertEquals(new AminoAcidSpecification(1, (String) null), aas(1, null));
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
        final TranscriptAminoAcids brafCanonical = baseSequenceVariantsCalculator.mTranscriptAminoAcidsMap.get(brafCanonicalTranscript.TransName);
        final int brafLength = brafCanonical.AminoAcids.length();
        assertTrue(aas(1, "M").applies(brafCanonical));
        assertTrue(aas(2, "A").applies(brafCanonical));
        assertTrue(aas(3, "A").applies(brafCanonical));
        assertFalse(aas(3, "L").applies(brafCanonical));
        assertTrue(aas(brafLength - 1, "H").applies(brafCanonical));
        assertTrue(aas(brafLength, "X").applies(brafCanonical));
        assertFalse(aas(brafLength + 1, "H").applies(brafCanonical));

        assertTrue(aas(4, null).applies(brafCanonical));
        assertTrue(aas(brafLength, null).applies(brafCanonical));
        assertFalse(aas(brafLength + 1, null).applies(brafCanonical));
    }
}
