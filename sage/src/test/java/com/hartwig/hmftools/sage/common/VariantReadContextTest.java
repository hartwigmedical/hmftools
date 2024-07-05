package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.Microhomology.findLeftHomologyShift;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.REF;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContextMatcher;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sage.candidate.AltRead;
import com.hartwig.hmftools.sage.candidate.RefContextConsumer;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class VariantReadContextTest
{
    private final static int TEST_FLANK_LENGTH = 5;

    @Test
    public void testSnvMnvVariantReadContext()
    {
        // test 1: most simple case with no repeats or homology
        SimpleVariant var = createSimpleVariant(50, "C", "A");

        String readBases = REF_BASES_200.substring(30, 50) + "A" + REF_BASES_200.substring(51, 70);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "40M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases, baseQuals);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(TEST_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(var, read, 20, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(43, readContext.AlignmentStart);
        assertEquals(57, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(7, readContext.VarIndex);
        assertEquals(9, readContext.CoreIndexEnd);
        assertEquals("CGCTG", readContext.refBases());
        assertEquals("GTCACCGATGTCTGT", readContext.readBases());
        assertEquals("15M", readContext.readCigar());
        assertEquals("GCT", readContext.trinucleotideStr());
        assertEquals(5, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("CGATG", readContext.coreStr());
        assertEquals("GTCAC", readContext.leftFlankStr());
        assertEquals("TCTGT", readContext.rightFlankStr());
        assertNull(readContext.Homology);
        assertNull(readContext.MaxRepeat);


        // test 2:
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

        readCigar = "15M3I10M3I15M";
        read = buildSamRecord(30, readCigar, readBases);

        readContext = builder.createContext(var, read, 23, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(44, readContext.AlignmentStart);
        assertEquals(55, readContext.AlignmentEnd);
        assertEquals(7, readContext.CoreIndexStart);
        assertEquals(9, readContext.VarIndex);
        assertEquals(11, readContext.CoreIndexEnd);
        assertEquals("CGCTG", readContext.refBases());
        assertEquals("TGGGCACCGATGTCGGGT", readContext.readBases());
        assertEquals("1M3I10M3I1M", readContext.readCigar());
        assertEquals(5, readContext.coreLength());
        assertEquals(7, readContext.leftFlankLength());
        assertEquals(6, readContext.rightFlankLength());
        assertEquals("CGATG", readContext.coreStr());
        assertEquals("TGGGCAC", readContext.leftFlankStr());
        assertEquals("TCGGGT", readContext.rightFlankStr());

        // test 3: with transitional repeats near the core

        // TGCGCGC TA CACACACT -> TGCGCGC GC CACACACT
        var = createSimpleVariant(143, "TA", "GC");

        readBases = REF_BASES_200.substring(120, 143) + "GC" + REF_BASES_200.substring(145, 170);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "50M";
        read = buildSamRecord(120, readCigar, readBases, baseQuals);

        readContext = builder.createContext(var, read, 23, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertNotNull(readContext.MaxRepeat);
        assertEquals("GC", readContext.MaxRepeat.Bases);
        assertEquals(4, readContext.MaxRepeat.Count);
        assertEquals(131, readContext.AlignmentStart);
        assertEquals(157, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(12, readContext.VarIndex);
        assertEquals(21, readContext.CoreIndexEnd);
        assertEquals("TGCGCGCTACACACACT", readContext.refBases());
        assertEquals("TGCGCGCGCCACACACT", readContext.coreStr());
        assertEquals("27M", readContext.readCigar());
        assertEquals(17, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
    }

    @Test
    public void testSnvComplex()
    {
        String refBases = REF_BASES_200.substring(0, 50) +
        //       50        60        70
        //       0        90123456789012345        4
        //                             T
                "CTTTTTTTTTATTTTTTTTTTTCTTTTTTTTGAGA" + REF_BASES_200.substring(100, 150);

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String readBases = refBases.substring(40, 72) + "T" + refBases.substring(73, 100);

        String readCigar = "60M";
        SAMRecord read = buildSamRecord(40, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SimpleVariant var = createSimpleVariant(72, "C", "T");
        VariantReadContext readContext = builder.createContext(var, read, 32, refSequence);

        assertTrue(readContext.isValid());
        assertEquals(50, readContext.AlignmentStart);
        assertEquals(91, readContext.AlignmentEnd);
        assertEquals(10, readContext.CoreIndexStart);
        assertEquals(22, readContext.VarIndex);
        assertEquals(31, readContext.CoreIndexEnd);
        assertEquals(12, readContext.variantRefIndex());
        assertEquals("ATTTTTTTTTTTCTTTTTTTTG", readContext.refBases());
        assertEquals("CTTTTTTTTTATTTTTTTTTTTTTTTTTTTTGAGAGAGGCAA", readContext.readBases());

        // test 1: a core match
        String readRefBases = refBases.substring(40, 100);
        SAMRecord refRead = buildSamRecord(40, readCigar, readRefBases);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);
        assertEquals(REF, matcher.determineReadMatch(refRead, 32));


        // test 2: repeats starting 1 base outside the basic core
        refBases = "AGGTCGATCGTAAAAAAATGCGTCGCGCGCGCGCGTAGGTCGATCG";
        //          100       110       120       130       140
        //          01234567890123456789012345678901234567890123456789
        //                              A
        refSequence = new RefSequence(100, refBases.getBytes());

        String altBase = "A";
        readBases = refBases.substring(0, 20) + "A" + refBases.substring(21);
        read = buildSamRecord(100, buildCigarString(readBases.length()), readBases);

        builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        var = createSimpleVariant(120, refBases.substring(20, 21), altBase);
        readContext = builder.createContext(var, read, 20, refSequence);

        assertTrue(readContext.isValid());
        assertEquals(100, readContext.AlignmentStart);
        assertEquals(145, readContext.AlignmentEnd);
        assertEquals(10, readContext.CoreIndexStart);
        assertEquals(20, readContext.VarIndex);
        assertEquals(35, readContext.CoreIndexEnd);
        assertEquals(2, readContext.AllRepeats.size());
        assertEquals("TAAAAAAATGAGTCGCGCGCGCGCGT", readContext.coreStr());
    }

    @Test
    public void testVariantRepeatContexts()
    {
        final String refBases =
            //             10        20        30        40        50
            //   0123456789012345678901234567890123456789012345678901
                "CGCAATATTCGGGTGGGAGTGACCCGATTTACCCGGTGCGTTCGTCACCGCTGTCT";

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        // test 1: a repeat starting at the variant
        String readBases = refBases.substring(0, 21) + "C" + refBases.substring(22, 51);

        String readCigar = buildCigarString(readBases.length());
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        SimpleVariant var = createSimpleVariant(21, "A", "C");

        VariantReadContextBuilder builder = new VariantReadContextBuilder(5);
        VariantReadContext readContext = builder.createContext(var, read, 21, refSequence);

        assertTrue(readContext.isValid());
        assertEquals(5, readContext.CoreIndexStart);
        assertNotNull(readContext.MaxRepeat);
        assertEquals("C", readContext.MaxRepeat.Bases);
        assertEquals(4, readContext.MaxRepeat.Count);
        assertEquals(7, readContext.VarIndex);
        assertEquals(11, readContext.CoreIndexEnd);

        // test 2 : repeats before and after but only the first is used
        var = createSimpleVariant(30, "A", "C");

        readBases = refBases.substring(0, 30) + var.alt() + refBases.substring(31, 51);

        readCigar = buildCigarString(readBases.length());
        read = buildSamRecord(1, readCigar, readBases);

        readContext = builder.createContext(var, read, 30, refSequence);

        assertTrue(readContext.isValid());
        assertEquals(5, readContext.CoreIndexStart);
        assertNotNull(readContext.MaxRepeat);
        assertEquals("C", readContext.MaxRepeat.Bases);
        assertEquals(4, readContext.MaxRepeat.Count);
        assertEquals(9, readContext.VarIndex);
        assertEquals(13, readContext.CoreIndexEnd);
    }

    @Test
    public void testDeleteVariantReadContext()
    {
        // firstly a simple delete will no repeats or homology
        SimpleVariant var = createSimpleVariant(50, "CTG", "C");

        String readBases = REF_BASES_200.substring(30, 51) + REF_BASES_200.substring(53, 70);
        String readCigar = "21M2D17M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(TEST_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(var, read, 20, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(44, readContext.AlignmentStart);
        assertEquals(60, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(6, readContext.VarIndex);
        assertEquals(9, readContext.CoreIndexEnd);
        assertEquals("GCTGTCT", readContext.refBases());
        assertEquals("TCACCGCTCTGTGAC", readContext.readBases());
        assertEquals("7M2D8M", readContext.readCigar());
        assertEquals(5, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("GCTCT", readContext.coreStr());
        assertEquals("TCACC", readContext.leftFlankStr());
        assertEquals("GTGAC", readContext.rightFlankStr());
        assertEquals(49, readContext.CorePositionStart);
        assertEquals(55, readContext.CorePositionEnd);

        // test 2: read context spans into the soft-clip bases on the left
        String refBases = REF_BASES_200.substring(57, 59);
        String altBase = refBases.substring(0, 1);
        var = createSimpleVariant(57, refBases, altBase);

        readBases = REF_BASES_200.substring(1, 58) + REF_BASES_200.substring(59, 79);

        // unclipped pos 1-49 = 49S, read index 0-48
        // aligned pos   50-57, read index 49-56
        // del pos       58
        // aligned pos   59-78, read index 57-76
        // so del occurs 49 + 8 = 57 bases into the read
        readCigar = "49S8M1D20M";
        read = buildSamRecord(50, readCigar, readBases);

        builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH); // full flank so it extends into the soft-clip
        readContext = builder.createContext(var, read, 56, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(50, readContext.AlignmentStart);
        assertEquals(70, readContext.AlignmentEnd);
        assertEquals(10, readContext.CoreIndexStart);
        assertEquals(11, readContext.VarIndex);
        assertEquals(13, readContext.CoreIndexEnd);
        assertEquals("4S8M1D12M", readContext.readCigar());
        assertEquals("GTAC", readContext.coreStr());
        assertEquals("ACCGCTGTCT", readContext.leftFlankStr());
        assertEquals("TCGGAAAAAA", readContext.rightFlankStr());

        // test 3: del spanning into the right soft-clip
        refBases = REF_BASES_200.substring(30, 33);
        altBase = refBases.substring(0, 1);
        var = createSimpleVariant(30, refBases, altBase);

        readBases = REF_BASES_200.substring(1, 31) + REF_BASES_200.substring(33, 53);
        readCigar = "30M2D5M15S";
        read = buildSamRecord(1, readCigar, readBases);

        builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH); // full flank so it extends into the soft-clip
        readContext = builder.createContext(var, read, 29, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(16, readContext.AlignmentStart);
        assertEquals(37, readContext.AlignmentEnd);
        assertEquals(10, readContext.CoreIndexStart);
        assertEquals(14, readContext.VarIndex);
        assertEquals(16, readContext.CoreIndexEnd);
        assertEquals("15M2D5M7S", readContext.readCigar());
    }

    @Test
    public void testInsertVariantReadContext()
    {
        SimpleVariant var = createSimpleVariant(50, "C", "CTG");

        String readBases = REF_BASES_200.substring(30, 51) + "TG" +  REF_BASES_200.substring(51, 70);
        String readCigar = "21M2I17M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(TEST_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(var, read, 20, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(44, readContext.AlignmentStart);
        assertEquals(60, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(6, readContext.VarIndex);
        assertEquals(13, readContext.CoreIndexEnd);
        // assertEquals(12, readContext.AltIndexUpper);
        assertEquals("GCTGTCT", readContext.refBases());
        assertEquals("TCACCGCTGTGTCTGTGAC", readContext.readBases());
        assertEquals("7M2I10M", readContext.readCigar());
        assertEquals(9, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("TG", readContext.homologyBases());
        assertEquals(3, readContext.Homology.Length);
        assertEquals("GCTGTGTCT", readContext.coreStr());
        assertEquals("TCACC", readContext.leftFlankStr());
        assertEquals("GTGAC", readContext.rightFlankStr());

        assertEquals(49, readContext.CorePositionStart);
        assertEquals(55, readContext.CorePositionEnd);


        // test 2: insert in left soft-clip with homology

        // a duplicated section leading up to the aligned bases
        // TTCGTCACCGCTGTCTGTGACTCGGAAAAAAAACTCC

        String insert = REF_BASES_200.substring(40, 55);
        readBases = REF_BASES_200.substring(15, 40) + insert +  REF_BASES_200.substring(40, 110);
        readCigar = "40S70M";
        read = buildSamRecord(40, readCigar, readBases);

        builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        AltRead altRead = RefContextConsumer.processSoftClip(
                40, 109, readBases, 40, 0, REF_SEQUENCE_200, true);

        var = createSimpleVariant(39, altRead.Ref, altRead.Alt);

        readContext = builder.createContext(var, read, 39, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(38, readContext.CorePositionStart);
        assertEquals(11, readContext.VarIndex);
        assertEquals(58, readContext.CorePositionEnd);
        assertEquals(16, readContext.Homology.Length);


        // test 3: insertion (a duplicated segment) in the right soft clip
        insert = REF_BASES_200.substring(85, 100);
        readBases = REF_BASES_200.substring(30, 100) + insert +  REF_BASES_200.substring(100, 125);
        readCigar = "70M40S";
        read = buildSamRecord(30, readCigar, readBases);

        builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        altRead = RefContextConsumer.processSoftClip(
                30, 99, readBases, 40, 70, REF_SEQUENCE_200, false);

        int readIndex = readBases.length() - 40 - 1;

        var = createSimpleVariant(99, altRead.Ref, altRead.Alt);

        int leftHomologyShift = findLeftHomologyShift(var, REF_SEQUENCE_200, readBases.getBytes(), readIndex);

        readIndex -= leftHomologyShift;
        int newPosition = var.Position - leftHomologyShift;
        String newRef = REF_SEQUENCE_200.positionBases(newPosition, newPosition);
        String newAlt = newRef + readBases.substring(readIndex + 1, readIndex + var.altLength());
        var = createSimpleVariant(newPosition, newRef, newAlt);

        readContext = builder.createContext(var, read, readIndex, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(80, readContext.CorePositionStart);
        assertEquals(14, readContext.VarIndex);
        assertEquals(116, readContext.CorePositionEnd);
    }
}
