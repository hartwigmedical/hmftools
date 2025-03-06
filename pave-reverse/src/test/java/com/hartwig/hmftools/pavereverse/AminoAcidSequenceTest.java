package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.List;

import org.junit.Test;

public class AminoAcidSequenceTest extends ReversePaveTestBase
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
        assertEquals(aa, aaSeq("AWR"));
        assertNotEquals(aa, aaSeq("AWRA"));
        assertNotEquals(aa, aaSeq("WR"));
    }

    @Test
    public void hashCodeTest()
    {
        List<AminoAcid> acids = List.of(aa("A"), aa("W"), aa("R"));
        AminoAcidSequence aa = new AminoAcidSequence(acids);
        assertEquals(aa.hashCode(), aaSeq("AWR").hashCode());
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
    public void replaceRangeTest()
    {
        assertEquals(aaSeq("SHECATSAREKIPPING"), aaSeq("THECATSAREKIPPING").replaceRange(1, 1, aaSeq("S")));
        assertEquals(aaSeq("THECARPAREKIPPING"), aaSeq("THECATSAREKIPPING").replaceRange(6, 7, aaSeq("RP")));
        assertEquals(aaSeq("THECATSAREKIPPINEES"), aaSeq("THECATSAREKIPPING").replaceRange(17, 17, aaSeq("EES")));
    }

    @Test
    public void includeStopCodonsWhenParsingNucleotidesTest()
    {
        assertEquals(aaSeq("X"), AminoAcidSequence.fromNucleotides("TGA"));
        assertEquals(aaSeq("AX"), AminoAcidSequence.fromNucleotides("GCTTAG"));
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

    @Test
    public void subsequenceUpToTest()
    {
        assertEquals(aaSeq(""), aaSeq("APCDEFGHI").subsequenceUpToInclusive(0));
        assertEquals(aaSeq("A"), aaSeq("APCDEFGHI").subsequenceUpToInclusive(1));
        assertEquals(aaSeq("AP"), aaSeq("APCDEFGHI").subsequenceUpToInclusive(2));
        assertEquals(aaSeq("APC"), aaSeq("APCDEFGHI").subsequenceUpToInclusive(3));
        assertEquals(aaSeq("APCDEFGH"), aaSeq("APCDEFGHI").subsequenceUpToInclusive(8));
        assertEquals(aaSeq("APCDEFGHI"), aaSeq("APCDEFGHI").subsequenceUpToInclusive(9));
    }

    @Test
    public void subsequenceAfterTest()
    {
        assertEquals(aaSeq(""), aaSeq("APCDEFGHI").subsequenceAfterExclusive(8));
        assertEquals(aaSeq("I"), aaSeq("APCDEFGHI").subsequenceAfterExclusive(7));
        assertEquals(aaSeq("HI"), aaSeq("APCDEFGHI").subsequenceAfterExclusive(6));
        assertEquals(aaSeq("CDEFGHI"), aaSeq("APCDEFGHI").subsequenceAfterExclusive(1));
        assertEquals(aaSeq("PCDEFGHI"), aaSeq("APCDEFGHI").subsequenceAfterExclusive(0));
    }

    @Test
    public void emptyTest()
    {
        assertEquals(0, AminoAcidSequence.empty().length());
    }

    @Test
    public void deleteTest()
    {
        assertEquals(aaSeq("HECATSAREKIPPING"), aaSeq("THECATSAREKIPPING").deleteRange(0, 1));
        assertEquals(aaSeq("THEKIPPING"), aaSeq("THECATSAREKIPPING").deleteRange(3, 10));
        assertEquals(aaSeq("G"), aaSeq("THECATSAREKIPPING").deleteRange(0, 16));
    }

    @Test
    public void duplicateRangeTest()
    {
        assertEquals(aaSeq("THECATSAREKIPPING"), aaSeq("THECATSAREKIPPING").duplicateRange(0, 0));
        assertEquals(aaSeq("THECATCATSAREKIPPING"), aaSeq("THECATSAREKIPPING").duplicateRange(3, 6));
        assertEquals(aaSeq("THECATSAREKIPPINGTHECATSAREKIPPING"), aaSeq("THECATSAREKIPPING").duplicateRange(0, 17));
    }

    @Test
    public void replaceTest()
    {
        assertEquals(aaSeq("SHECATSAREKIPPING"), aaSeq("THECATSAREKIPPING").replace(1, aa("S")));
        assertEquals(aaSeq("THECADSAREKIPPING"), aaSeq("THECATSAREKIPPING").replace(6, aa("D")));
        assertEquals(aaSeq("THECATSAREKIPPINE"), aaSeq("THECATSAREKIPPING").replace(17, aa("E")));
    }

    @Test
    public void insertTest()
    {
        assertEquals(aaSeq("ALLTHECATSAREKIPPING"), aaSeq("THECATSAREKIPPING").insert(0, aaSeq("ALL")));
        assertEquals(aaSeq("THEGIANTCATSAREKIPPING"), aaSeq("THECATSAREKIPPING").insert(3, aaSeq("GIANT")));
        assertEquals(aaSeq("THECATSAREKIPPINGHERE"), aaSeq("THECATSAREKIPPING").insert(17, aaSeq("HERE")));
    }

    @Test
    public void appendTest()
    {
        assertEquals(aaSeq("APCDEFGHI"), aaSeq("APCDE").append(aaSeq("FGHI")));
    }
}
