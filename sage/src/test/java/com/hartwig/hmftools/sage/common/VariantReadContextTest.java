package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.common.VariantReadContextBuilder.findHomology;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class VariantReadContextTest
{
    @Test
    public void testSnvMnvVariantReadContext()
    {
        SimpleVariant var = createSimpleVariant(50, "C", "A");

        String readBases = REF_BASES_200.substring(30, 50) + "A" + REF_BASES_200.substring(51, 70);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "40M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases, baseQuals);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(5);

        VariantReadContext readContext = builder.createMnvContext(var, read, 20, REF_SEQUENCE_200);

        assertEquals(43, readContext.AlignmentStart);
        assertEquals(57, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(7, readContext.VarReadIndex);
        assertEquals(9, readContext.CoreIndexEnd);
        assertEquals("GTCACCGCTGTCTGT", readContext.refBases());
        assertEquals("GTCACCGATGTCTGT", readContext.readBases());
        assertEquals("15M", readContext.readCigar());
        assertEquals("GCT", readContext.trinucleotideStr());
        assertEquals(5, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("CGATG", readContext.coreStr());
        assertEquals("GTCAC", readContext.leftFlankStr());
        assertEquals("TCTGT", readContext.rightFlankStr());

        // Read bases and cigar:
        // 15M index 0-14, pos 30-44
        // 3I index 15-17, from pos 44
        // 10M index 18-27, pos 45-54, with SNV at index 23
        // 3I index 28-30, from pos 54
        // 15M index 31-45, pos 55-69
        // Initially flank indices are 16-30 of the read bases, but since 16 falls into the insert, it will push back to 14 and 1M at 44
        // and likewise since 30 falls into the last insert it will push back to 31 and pos 55

        readBases = REF_BASES_200.substring(30, 45) + "GGG" + REF_BASES_200.substring(45, 50) + "A"
                + REF_BASES_200.substring(51, 55) + "GGG" + REF_BASES_200.substring(55, 70);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "15M3I10M3I15M";
        read = buildSamRecord(30, readCigar, readBases, baseQuals);

        readContext = builder.createMnvContext(var, read, 23, REF_SEQUENCE_200);

        assertEquals(44, readContext.AlignmentStart);
        assertEquals(55, readContext.AlignmentEnd);
        assertEquals(7, readContext.CoreIndexStart);
        assertEquals(9, readContext.VarReadIndex);
        assertEquals(11, readContext.CoreIndexEnd);
        assertEquals("TCACCGCTGTCT", readContext.refBases());
        assertEquals("TGGGCACCGATGTCGGGT", readContext.readBases());
        assertEquals("1M2I10M3I1M", readContext.readCigar());
        assertEquals(5, readContext.coreLength());
        assertEquals(7, readContext.leftFlankLength());
        assertEquals(6, readContext.rightFlankLength());
        assertEquals("CGATG", readContext.coreStr());
        assertEquals("TGGGCAC", readContext.leftFlankStr());
        assertEquals("TCGGGT", readContext.rightFlankStr());

        // with transitional repeats near the core

        // TGCGCGC TA CACACACT -> TGCGCGC GC CACACACT
        var = createSimpleVariant(143, "TA", "GC");

        readBases = REF_BASES_200.substring(120, 143) + "GC" + REF_BASES_200.substring(145, 170);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "50M";
        read = buildSamRecord(120, readCigar, readBases, baseQuals);

        readContext = builder.createMnvContext(var, read, 23, REF_SEQUENCE_200);

        assertEquals(131, readContext.AlignmentStart);
        assertEquals(157, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(12, readContext.VarReadIndex);
        assertEquals(21, readContext.CoreIndexEnd);
        assertEquals("CACAGTGCGCGCTACACACACTGGCCT", readContext.refBases());
        assertEquals("CACAGTGCGCGCGCCACACACTGGCCT", readContext.readBases());
        assertEquals("27M", readContext.readCigar());
        assertEquals(17, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("TGCGCGCGCCACACACT", readContext.coreStr());
    }

    @Test
    public void testDeleteVariantReadContext()
    {
        // firstly a simple delete will no repeats or homology
        SimpleVariant var = createSimpleVariant(50, "CTG", "C");

        String readBases = REF_BASES_200.substring(30, 51) + REF_BASES_200.substring(53, 70);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "21M2D17M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases, baseQuals);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(5);

        VariantReadContext readContext = builder.createIndelContext(var, read, 20, REF_SEQUENCE_200);

        assertEquals(44, readContext.AlignmentStart);
        assertEquals(60, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(6, readContext.VarReadIndex);
        assertEquals(9, readContext.CoreIndexEnd);
        assertEquals("TCACCGCTGTCTGTGAC", readContext.refBases());
        assertEquals("TCACCGCTCTGTGAC", readContext.readBases());
        assertEquals("7M2D8M", readContext.readCigar());
        assertEquals(5, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("GCTCT", readContext.coreStr());
        assertEquals("TCACC", readContext.leftFlankStr());
        assertEquals("GTGAC", readContext.rightFlankStr());
    }

    @Test
    public void testHomology()
    {
        SimpleVariant var = createSimpleVariant(26, "A", "AT");

        String readBases = REF_BASES_200.substring(10, 26) + "AT" + REF_BASES_200.substring(27, 40);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "16M1I13M";
        SAMRecord read = buildSamRecord(10, readCigar, readBases, baseQuals);

        Microhomology homology = findHomology(var, read, 16, REF_SEQUENCE_200);

        assertNotNull(homology);
        assertEquals("T", homology.Bases);
        assertEquals(4, homology.Length);

        var = createSimpleVariant(26, "A", "ATT");

        readBases = REF_BASES_200.substring(10, 26) + "ATT" + REF_BASES_200.substring(27, 40);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "16M2I13M";
        read = buildSamRecord(10, readCigar, readBases, baseQuals);

        homology = findHomology(var, read, 16, REF_SEQUENCE_200);

        assertNotNull(homology);
        assertEquals("TT", homology.Bases);
        assertEquals(4, homology.Length);

        // delete
        var = createSimpleVariant(64, "GAAA", "G");

        readBases = REF_BASES_200.substring(50, 65) + REF_BASES_200.substring(69, 80);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "15M3D11M";
        read = buildSamRecord(50, readCigar, readBases, baseQuals);

        homology = findHomology(var, read, 15, REF_SEQUENCE_200);

        assertNotNull(homology);
        assertEquals("AAA", homology.Bases);
        assertEquals(3, homology.Length);
    }
}
