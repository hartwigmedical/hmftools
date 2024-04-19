package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.RepeatInfo.findMultiBaseRepeat;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantReadContextBuilder.findHomology;
import static com.hartwig.hmftools.sage.common.VariantUtils.TEST_LEFT_CORE;
import static com.hartwig.hmftools.sage.common.VariantUtils.TEST_RIGHT_CORE;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class VariantReadContextTest
{
    private final static int TEST_FLANK_LENGTH = 5;

    @Test
    public void testSnvMnvVariantReadContext()
    {
        SimpleVariant var = createSimpleVariant(50, "C", "A");

        String readBases = REF_BASES_200.substring(30, 50) + "A" + REF_BASES_200.substring(51, 70);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "40M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases, baseQuals);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(TEST_FLANK_LENGTH);

        VariantReadContext readContext = builder.createSnvMnvContext(var, read, 20, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(43, readContext.AlignmentStart);
        assertEquals(57, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(7, readContext.VarReadIndex);
        assertEquals(9, readContext.CoreIndexEnd);
        assertEquals(7, readContext.AltIndexLower);
        assertEquals(7, readContext.AltIndexUpper);
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

        readContext = builder.createSnvMnvContext(var, read, 23, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(44, readContext.AlignmentStart);
        assertEquals(55, readContext.AlignmentEnd);
        assertEquals(7, readContext.CoreIndexStart);
        assertEquals(9, readContext.VarReadIndex);
        assertEquals(11, readContext.CoreIndexEnd);
        assertEquals(9, readContext.AltIndexLower);
        assertEquals(9, readContext.AltIndexUpper);
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

        readContext = builder.createSnvMnvContext(var, read, 23, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(131, readContext.AlignmentStart);
        assertEquals(157, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(12, readContext.VarReadIndex);
        assertEquals(21, readContext.CoreIndexEnd);
        assertEquals(12, readContext.AltIndexLower);
        assertEquals(13, readContext.AltIndexUpper);
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

        VariantReadContextBuilder builder = new VariantReadContextBuilder(TEST_FLANK_LENGTH);

        VariantReadContext readContext = builder.createIndelContext(var, read, 20, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(44, readContext.AlignmentStart);
        assertEquals(60, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(6, readContext.VarReadIndex);
        assertEquals(9, readContext.CoreIndexEnd);
        assertEquals(6, readContext.AltIndexLower);
        assertEquals(7, readContext.AltIndexUpper);
        assertEquals("TCACCGCTGTCTGTGAC", readContext.refBases());
        assertEquals("TCACCGCTCTGTGAC", readContext.readBases());
        assertEquals("7M2D8M", readContext.readCigar());
        assertEquals(5, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("GCTCT", readContext.coreStr());
        assertEquals("TCACC", readContext.leftFlankStr());
        assertEquals("GTGAC", readContext.rightFlankStr());
        assertEquals(49, readContext.corePositionStart());
        assertEquals(55, readContext.corePositionEnd());

        // CLEAN-UP: add more scenarios
    }

    @Test
    public void testInsertVariantReadContext()
    {
        SimpleVariant var = createSimpleVariant(50, "C", "CTG");

        String readBases = REF_BASES_200.substring(30, 51) + "TG" +  REF_BASES_200.substring(51, 70);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "21M2I17M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases, baseQuals);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(TEST_FLANK_LENGTH);

        VariantReadContext readContext = builder.createIndelContext(var, read, 20, REF_SEQUENCE_200);

        assertTrue(readContext.isValid());
        assertEquals(44, readContext.AlignmentStart);
        assertEquals(60, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(6, readContext.VarReadIndex);
        assertEquals(13, readContext.CoreIndexEnd);
        assertEquals(6, readContext.AltIndexLower);
        assertEquals(9, readContext.AltIndexUpper);
        assertEquals("TCACCGCTGTCTGTGAC", readContext.refBases());
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

        assertEquals(49, readContext.corePositionStart());
        assertEquals(55, readContext.corePositionEnd());

        // CLEAN-UP: add more scenarios
    }

    /* old scenarios:

        @Test
        public void testSimpleSnvHas5BaseCore()
        {
            String refSequence = "GATCATCTAGG";
            String readSequence = "GATCACCTAGG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("11M", readSequence);
            ReadContext victim = this.victim.createSNVContext(1005, 5, record, refBases);
            assertEquals("CACCT", victim.coreString());
        }

        @Test
        public void testSimpleInsert()
        {
            String refSequence = "GATCATCTAGG";
            String readSequence = "GAGGCTCATCTAGG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("2M3I9M", readSequence);
            ReadContext victim = this.victim.createInsertContext("AGGC", 1000, 1, record.getReadBases(), refBases);
            assertEquals("GAGGCTC", victim.coreString());
        }

        @Test
        public void testInsertInRepeat()
        {
            String refSequence = "TGAAAAAAAATCT";
            String readSequence = "TGAAAAAAAAATCT";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("2M1I11M", readSequence);
            ReadContext victim = this.victim.createInsertContext("GA", 1000, 1, record.getReadBases(), refBases);
            assertEquals("TGAAAAAAAAAT", victim.coreString());
        }

        @Test
        public void testInsertAtHomology()
        {
            String refSequence = "GATCATCTG";
            String readSequence = "GATCATCATCTG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("1M3I8M", readSequence);
            ReadContext victim = this.victim.createInsertContext("ATCA", 1000, 1, record.getReadBases(), refBases);
            assertEquals("GATCATCATCT", victim.coreString());
        }

        @Test
        public void testInsertAtHomologyRepeat()
        {
            String refSequence = "GATCATCATCTG";
            String readSequence = "GATCATCATCATCTG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("1M3I10M", readSequence);
            ReadContext victim = this.victim.createInsertContext("ATCA", 1000, 1, record.getReadBases(), refBases);
            assertEquals("GATCATCATCATCT", victim.coreString());
        }

        @Test
        public void testInsertAtHomologyWithAdditionalBases()
        {
            String refSequence = "ATGCGATCTTCC";
            String readSequence = "ATGCGATCAATCTTCC";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("5M4I7M", readSequence);
            ReadContext victim = this.victim.createInsertContext("GATCA", 1000, 4, record.getReadBases(), refBases);
            assertEquals("GCGATCAAT", victim.coreString());
        }

        private static final ReadContextFactory READ_CONTEXT_FACTORY = new ReadContextFactory(DEFAULT_FLANK_LENGTH);

        @Test
        public void testSimpleDelete1()
        {
            // variant: pos 1025 AT>A
            String refBases =  "ATCTCTCAATGTTGACGGACAGCCTATTTTTGCCAATATCACACTGCCAGGT";
            String readBases = "ATCTCTCAATGTTGACGGACAGCCTATTTTGCCAATATCACACTGCCAGGT";
            RefSequence refSequence = new RefSequence(1000, refBases.getBytes());

            ReadContext readContext = READ_CONTEXT_FACTORY.createDelContext("AT", 1025, 25, readBases.getBytes(), refSequence);

            assertFalse(readContext.hasIncompleteCore());
            assertEquals("CTATTTTGC", readContext.coreString());
            assertEquals("GACGGACAGC", readContext.leftFlankString());
            assertEquals("CAATATCACA", readContext.rightFlankString());
        }

        @Test
        public void testDeleteAtHomology()
        {
            String refSequence = "GATCGGATCGCTT";
            String readSequence = "GATCGCTT";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("1M5D7M", readSequence);
            ReadContext victim = this.victim.createDelContext("GATCGG", 1000, 0, record.getReadBases(), refBases);
            assertEquals("GATCGC", victim.coreString());
        }

        @Test
        public void testDeleteAtHomologyRepeat()
        {
            String refSequence = "GATCACCATCTG";
            String readSequence = "GATCATCTG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("1M3D8M", readSequence);
            ReadContext victim = this.victim.createDelContext("ATCA", 1000, 1, record.getReadBases(), refBases);
            assertEquals("GATCATCT", victim.coreString());
        }

        @Test
        public void testDeleteAtRepeatInRef()
        {
            String refSequence = "GATCATCATCTG";
            String readSequence = "GATCATCTG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("1M3D8M", readSequence);
            ReadContext victim = this.victim.createDelContext("ATCA", 1000, 1, record.getReadBases(), refBases);
            assertEquals("GATCATCTG", victim.coreString());
        }

        @Test
        public void testDeleteOneBase()
        {
            String refSequence = "GATCATCTAGG";
            String readSequence = "GTCATCTAGG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("1M1D9M", readSequence);
            ReadContext victim = this.victim.createDelContext("GA", 1000, 0, record.getReadBases(), refBases);
            assertEquals("GTC", victim.coreString());
        }

        @Test
        public void testDeleteTwoBase()
        {
            String refSequence = "GATCATCTAGG";
            String readSequence = "GCATCTAGG";
            RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

            SAMRecord record = buildSamRecord("1M2D8M", readSequence);
            ReadContext victim = this.victim.createDelContext("GAT", 1000, 0, record.getReadBases(), refBases);
            assertEquals("GCA", victim.coreString());
        }


     */

    @Test
    public void testRepeats()
    {
        //              0123456789
        String bases = "AAACCTTTTT";

        // first check limits
        RepeatInfo repeatInfo = findMultiBaseRepeat(bases.getBytes(), 6, 1, 4);
        assertNotNull(repeatInfo);
        assertEquals("T", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);

        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 7, 1, 4);
        assertNull(repeatInfo);

        //       01234567890
        bases = "AACGGTACGGT";

        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 1, 5, 2);
        assertNotNull(repeatInfo);
        assertEquals("ACGGT", repeatInfo.Bases);
        assertEquals(2, repeatInfo.Count);
        assertEquals(10, repeatInfo.length());
        assertEquals(5, repeatInfo.repeatLength());

        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 2, 5, 2);
        assertNull(repeatInfo);
    }

    @Test
    public void testHomology()
    {
        SimpleVariant var = createSimpleVariant(26, "A", "AT");

        String readBases = REF_BASES_200.substring(10, 26) + "AT" + REF_BASES_200.substring(27, 40);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "16M1I13M";
        SAMRecord read = buildSamRecord(10, readCigar, readBases, baseQuals);

        Microhomology homology = findHomology(var, read, 16);

        assertNotNull(homology);
        assertEquals("T", homology.Bases);
        assertEquals(4, homology.Length);

        var = createSimpleVariant(26, "A", "ATT");

        readBases = REF_BASES_200.substring(10, 26) + "ATT" + REF_BASES_200.substring(27, 40);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "16M2I13M";
        read = buildSamRecord(10, readCigar, readBases, baseQuals);

        homology = findHomology(var, read, 16);

        assertNotNull(homology);
        assertEquals("TT", homology.Bases);
        assertEquals(4, homology.Length);

        // delete
        var = createSimpleVariant(64, "GAAA", "G");

        readBases = REF_BASES_200.substring(50, 65) + REF_BASES_200.substring(69, 80);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "15M3D11M";
        read = buildSamRecord(50, readCigar, readBases, baseQuals);

        homology = findHomology(var, read, 15);

        assertNotNull(homology);
        assertEquals("AAA", homology.Bases);
        assertEquals(3, homology.Length);

        // checks read bases not ref bases to determine homology
        var = createSimpleVariant(26, "A", "AAAA");
        readBases = REF_BASES_200.substring(10, 26) + "AAAAAAAAAAAAAAAAAAAAATG";
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "16M3I20M";
        read = buildSamRecord(10, readCigar, readBases, baseQuals);

        homology = findHomology(var, read, 16);

        assertNotNull(homology);
        assertEquals("AAA", homology.Bases);
        assertEquals(17, homology.Length);

        // multiple copies and then finishes with a partial copy
        // eg G(AACTC)AACTCAACTCAACCCTTT -> GAACTCAACTCAA(CTCAA)CCCTTT

        var = createSimpleVariant(26, "G", "GAACTC");
        readBases = REF_BASES_200.substring(10, 26) + "GAACTCAACTCAACTCAACCCTTT";
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "16M5I18M";
        read = buildSamRecord(10, readCigar, readBases, baseQuals);

        homology = findHomology(var, read, 16);

        assertNotNull(homology);
        assertEquals("AACTC", homology.Bases);
        assertEquals(13, homology.Length);
    }
}
