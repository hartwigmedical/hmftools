package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.List;

import org.junit.Test;

public class AminoAcidSequenceTest extends TransvalTestBase
{
    @Test
    public void sequenceTest()
    {
        List<AminoAcid> acids = List.of(aa("A"), aa("W"), aa("R"));
        assertEquals("AWR", new AminoAcidSequence(acids).sequence());
    }

    @Test
    public void parseTest()
    {
        assertEquals("A", aaSeq("A").sequence());
        assertEquals("A", aaSeq("Ala").sequence());
        assertEquals("R", aaSeq("Arg").sequence());
        assertEquals("RA", aaSeq("ArgA").sequence());
        assertEquals("YWV", aaSeq("YWV").sequence());
        assertEquals("YWV", aaSeq("TyrWV").sequence());
        assertEquals("YWV", aaSeq("TyrWVal").sequence());
        assertEquals("YWV", aaSeq("TyrTrpVal").sequence());
    }

    @Test
    public void equalsTest()
    {
        List<AminoAcid> acids = List.of(aa("A"), aa("W"), aa("R"));
        AminoAcidSequence aa = new AminoAcidSequence(acids);
        assertEquals(aa,  aaSeq("AWR"));
        assertNotEquals(aa,  aaSeq("AWRA"));
        assertNotEquals(aa,  aaSeq("WR"));
    }

    @Test
    public void hashCodeTest()
    {
        List<AminoAcid> acids = List.of(aa("A"), aa("W"), aa("R"));
        AminoAcidSequence aa = new AminoAcidSequence(acids);
        assertEquals(aa.hashCode(),  aaSeq("AWR").hashCode());
    }

    @Test
    public void toStringTest()
    {
        assertEquals("AWQ", aaSeq("AWQ").toString());
    }

    @Test
    public void fromNucleotidesTest()
    {
        assertEquals(aaSeq("A"), AminoAcidSequence.fromNucleotides("GCT"));
        assertEquals(aaSeq("AA"), AminoAcidSequence.fromNucleotides("GCTGCT"));
        assertEquals(aaSeq("AA"), AminoAcidSequence.fromNucleotides("GCTGCA"));
        assertEquals(aaSeq("SA"), AminoAcidSequence.fromNucleotides("AGTGCA"));
        assertEquals(aaSeq("SAT"), AminoAcidSequence.fromNucleotides("AGTGCAACC"));
    }

    @Test
    public void symbolAtTest()
    {
        assertEquals("A", aaSeq("A").symbolAt(0));
    }

    @Test
    public void lengthTest()
    {
        assertEquals(1, aaSeq("A").length());
        assertEquals(3, aaSeq("SAT").length());
    }
}
