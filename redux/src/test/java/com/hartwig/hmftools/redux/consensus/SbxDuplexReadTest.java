package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_YC_TAG;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.stripDuplexIndels;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.common.FragmentCoords;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SbxDuplexReadTest
{
    private static final String MISMATCH_QUAL = String.valueOf(phredToFastq(SBX_DUPLEX_MISMATCH_QUAL));
    private static final String NON_ZERO_QUAL = String.valueOf(phredToFastq(RAW_DUPLEX_QUAL));

    public SbxDuplexReadTest()
    {
        SbxRoutines.SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL = false;
    }

    @After
    public void resetSequencingType()
    {
        SbxRoutines.SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL = true;
    }

    @Test
    public void testHomopolymerLowBaseQualDirection()
    {
        //                  01234567890123456789
        String readBases = "ACGTACGTAAAAACGTACGT";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, 100, readBases, "20M", false, false, null);

        byte[] baseQualities = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        read.setBaseQualities(baseQualities);

        assertNull(SbxBamUtils.isHomopolymerLowBaseQualAtStart(read));

        read.getBaseQualities()[8] = 0;
        assertTrue(SbxBamUtils.isHomopolymerLowBaseQualAtStart(read));

        read.getBaseQualities()[8] = RAW_DUPLEX_QUAL;
        read.getBaseQualities()[11] = 0;
        read.getBaseQualities()[12] = 0;
        assertFalse(SbxBamUtils.isHomopolymerLowBaseQualAtStart(read));

        // unclear on a single base
        read.getBaseQualities()[11] = RAW_DUPLEX_QUAL;
        read.getBaseQualities()[12] = RAW_DUPLEX_QUAL;
        read.getBaseQualities()[15] = 0;
        assertNull(SbxBamUtils.isHomopolymerLowBaseQualAtStart(read));
    }

    @Test
    public void testStripDuplexIndels()
    {
        int alignmentStart = 10;

        // test invalid: no stripping for a SNV

        //                012567
        String readStr = "CTCACC";
        String cigar = "6M";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        byte[] baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[3] = 0;
        read.setBaseQualities(baseQuals);

        String ycTagStr = "0-3Z2-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 1);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("6M", read.getCigarString());
        assertEquals("CTCACC", read.getReadString());
        checkBaseQuals(read, List.of(3));
        assertEquals(1, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 0: 1 low-qual base not in a repeat

        //         01234567
        readStr = "CTCTACC";
        cigar = "3M1I3M";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[3] = 0;
        read.setBaseQualities(baseQuals);

        ycTagStr = "0-3Z3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 1);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("6M", read.getCigarString());
        assertEquals("CTCACC", read.getReadString());
        checkBaseQuals(read, Collections.emptyList());
        assertEquals(0, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 2: as before in 2-base repeat

        //         012345678
        readStr = "CTCTTTACC";
        cigar = "3M2I4M";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[5] = 0;
        read.setBaseQualities(baseQuals);

        ycTagStr = "0-5Z3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("3M1I4M", read.getCigarString());
        assertEquals("CTCTTACC", read.getReadString());
        checkBaseQuals(read, List.of(3, 4));
        assertEquals(1, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 3: 2 low-qual bases and 2-base repeat

        //         0123456789
        readStr = "CTCTTTTACC";
        cigar = "3M2I5M";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[5] = 0;
        baseQuals[6] = 0;
        read.setBaseQualities(baseQuals);

        ycTagStr = "0-5ZZ3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("8M", read.getCigarString());
        assertEquals("CTCTTACC", read.getReadString());
        checkBaseQuals(read, List.of(3, 4));
        assertEquals(0, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 4: 2 low-qual bases and 2-base repeat

        //         012     345     67890
        readStr = "CTC" + "TTT" + "TTACC";
        cigar = "3M3I5M";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[5] = 0;
        baseQuals[6] = 0;
        read.setBaseQualities(baseQuals);

        ycTagStr = "0-6ZZ3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 3);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("3M1I5M", read.getCigarString());
        assertEquals("CTCTTTACC", read.getReadString());
        checkBaseQuals(read, List.of(3, 5));
        assertEquals(1, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 5: 2 low-qual bases and 2-base repeat with A not part of the repeat

        //         012     3456     78901
        readStr = "CTC" + "ATTT" + "TTACC";
        cigar = "3M4I5M";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[6] = 0;
        baseQuals[7] = 0;
        read.setBaseQualities(baseQuals);

        ycTagStr = "0-7ZZ3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 4);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("3M2I5M", read.getCigarString());
        assertEquals("CTCATTTACC", read.getReadString());
        checkBaseQuals(read, List.of(4, 6));
        assertEquals(2, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
    }

    @Test
    public void testStripDuplexIndelsComplex()
    {
        int alignmentStart = 10;

        // test 1: Non-homopolymer repeat

        //                012     3456     78901
        String readStr = "CTC" + "ATAT" + "ATACC";
        String cigar = "3M4I5M";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        byte[] baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[7] = 0;
        baseQuals[8] = 0;
        read.setBaseQualities(baseQuals);

        String ycTagStr = "0-7ZZ3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 4);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("3M2I5M", read.getCigarString());
        assertEquals("CTCATATACC", read.getReadString());
        checkBaseQuals(read, List.of(3, 6));
        assertEquals(2, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 2: requiring left alignment after trimming
        //         012     3456     78     901
        readStr = "CTC" + "AAAA" + "AG" + "CTC";
        cigar = "7M2I3M";

        // M MMMM II M    M MMMM I M    M I MMMM M
        // C AAAA AG C -> C AAAA A C -> C A AAAA C
        // 9 9999 90 9    9 0999 0 9    9 0 9990 9

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[8] = 0;
        read.setBaseQualities(baseQuals);

        ycTagStr = "0-8Z3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        // read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("3M1I7M", read.getCigarString());
        assertEquals("CTCAAAAACTC", read.getReadString());
        checkBaseQuals(read, Collections.emptyList());
        assertEquals(1, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 3: low-qual base is right aligned
        // MMMMM I MMMMMMMM
        // CTTTA T TTTTATTA
        // 99990 9 99909999
        //         01234     5     6789012
        readStr = "CTTTA" + "T" + "TTTTACC";
        cigar = "5M1I7M";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[4] = 0;
        baseQuals[9] = 0;
        read.setBaseQualities(baseQuals);

        ycTagStr = "0-4Z4Z3-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        // read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 10);

        stripDuplexIndels(read);

        assertEquals("12M", read.getCigarString());
        assertEquals("CTTTATTTTACC", read.getReadString());
        checkBaseQuals(read, List.of(4, 5, 8));
        assertEquals(1, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // test 4: multiple stripped indels need to offset each other
        //         0123     456     7690     12     3456     7     8901
        readStr = "AACC" + "TTT" + "AACC" + "TT" + "AACC" + "T" + "AACC";
        cigar = "4M3I4M2I4M1I4M";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        baseQuals = buildBaseQuals(readStr.length(), RAW_DUPLEX_QUAL);
        baseQuals[4] = 0;
        baseQuals[5] = 0;
        baseQuals[6] = 0;

        baseQuals[11] = 0;
        baseQuals[12] = 0;

        baseQuals[17] = 0;

        read.setBaseQualities(baseQuals);

        ycTagStr = "0-4ZZZ4ZZ4Z4-0";

        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 6);

        stripDuplexIndels(read);

        assertEquals("16M", read.getCigarString());
        assertEquals("AACCAACCAACCAACC", read.getReadString());
        checkBaseQuals(read, Collections.emptyList());
        assertEquals(0, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
   }

    private static void checkBaseQuals(final SAMRecord read, final List<Integer> lowQualBases)
    {
        for(int i = 0; i < read.getBaseQualities().length; ++i)
        {
            if(lowQualBases.contains(i))
            {
                if(read.getBaseQualities()[i] != RAW_DUPLEX_MISMATCH_QUAL && read.getBaseQualities()[i] != SBX_DUPLEX_MISMATCH_QUAL)
                    assertFalse(true);
            }
            else
            {
                assertEquals(RAW_DUPLEX_QUAL, read.getBaseQualities()[i]);
            }
        }
    }

    @Test
    public void testStripDuplexIndelsForwardRead()
    {
        int alignmentStart = 25;

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(100);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        int mapq = 10;
        int nm = 2;
        int alignmentScore = 0;
        String readStr = "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        String cigar = "5M2I7M";
        String qualStr = NON_ZERO_QUAL.repeat(7) + MISMATCH_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);
        String ycTagStr = "0-7ZZ5-0";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(mapq);
        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, nm);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignmentScore);

        stripDuplexIndels(read);

        int alignmentScoreDiff = BWA_GAP_OPEN_PENALTY + 2 * BWA_GAP_EXTEND_PENALTY;

        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(5) + "CT" + "A".repeat(5);
        String expectedCigar = "12M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(5) + MISMATCH_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);

        SAMRecord expectedRead = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, expectedReadStr, expectedCigar, false, false, null);

        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(mapq);
        expectedRead.setAttribute(SBX_YC_TAG, ycTagStr);
        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, expectedNm);
        expectedRead.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, expectedAlignmentScore);

        assertEquals(expectedRead, read);
    }

    @Test
    public void testStripDuplexIndelsReverseRead()
    {
        int alignmentStart = 25;

        int mapq = 10;
        int nm = 2;
        int alignmentScore = 0;
        String readStr = "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        String cigar = "5M2I7M";  // inserts are left aligned.
        String qualStr = NON_ZERO_QUAL.repeat(5) + MISMATCH_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(7);
        String ycTagStr = "0-7ZZ5-0";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, true, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(mapq);
        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, nm);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignmentScore);

        stripDuplexIndels(read);

        int alignmentScoreDiff = BWA_GAP_OPEN_PENALTY + 2 * BWA_GAP_EXTEND_PENALTY;

        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(5) + "CT" + "A".repeat(5);
        String expectedCigar = "12M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(5) + MISMATCH_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);

        SAMRecord expectedRead = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, expectedReadStr, expectedCigar, true, false, null);

        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(mapq);
        expectedRead.setAttribute(SBX_YC_TAG, ycTagStr);
        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, expectedNm);
        expectedRead.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, expectedAlignmentScore);

        assertEquals(expectedRead, read);
    }

    @Test
    public void testStripDuplexIndelsInWithSupplementary()
    {
        // use supplementary alignment data to strip indels in soft-clip region
        int alignmentStart = 25;

        String alignedBases = REF_BASES.substring(1, 21);

        // index:           20             30
        //                  0123456789     0     1234567890
        String suppBases = "ACGTACGTAA" + "G" + "ACGTACGTAA";
        String cigar = "20M21S";
        String suppCigar = "20S10M1I10M";
        int suppPosStart = 25;
        String readBases = alignedBases + suppBases;

        SupplementaryReadData suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_POS_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        String ycTagStr = "20-10Z13-0";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);
        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.getBaseQualities()[30] = RAW_DUPLEX_MISMATCH_QUAL;

        stripDuplexIndels(read);

        assertEquals("20M20S", read.getCigarString());
        assertEquals(40, read.getReadBases().length);

        /*
        // with supplementary having a different aligned and soft-clip length, requiring adjustment
        cigar = "20M21S";
        suppCigar = "18S12M2I12M";

        suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_POS_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);
        read.setAttribute(SBX_YC_TAG, ycTagStr);

        stripDuplexIndels(read);

        assertEquals("20M23S", read.getCigarString());
        assertEquals(43, read.getReadBases().length);
        */

        // reversed supp data
        // index:    20                  30
        //           0123     4     5678901234
        suppBases = "ACGT" + "G" + "ACGTACGTAA";    // REF_BASES.substring(30, 54); // "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        readBases = alignedBases + suppBases;

        cigar = "20M15S";
        suppCigar = "10M1I4M";
        suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_NEG_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);

        ycTagStr = "0-24Z10-0";
        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.getBaseQualities()[24] = RAW_DUPLEX_MISMATCH_QUAL;

        stripDuplexIndels(read);

        assertEquals("20M14S", read.getCigarString());
        assertEquals(34, read.getReadBases().length);

        // repeat with soft-clip on the other side
        cigar = "21S20M";
        suppCigar = "10M1I10M20S";
        suppBases = "ACGTACGTAA" + "G" + "ACGTACGTAA";
        readBases = suppBases + alignedBases;

        suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_POS_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        ycTagStr = "0-10Z30-0";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);
        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.getBaseQualities()[10] = RAW_DUPLEX_MISMATCH_QUAL;

        stripDuplexIndels(read);

        assertEquals("20S20M", read.getCigarString());
        assertEquals(40, read.getReadBases().length);
    }

    private static RefGenomeInterface getRefGenome(final String refBases)
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, refBases);
        mockRefGenome.ChromosomeLengths.put(CHR_1, refBases.length());
        return new CachedRefGenome(mockRefGenome);
    }

    // UNCLEAR of the value of these tests
    private static final String DUPLEX_QUAL_STR = String.valueOf(phredToFastq(RAW_DUPLEX_QUAL));

    public static final String TEST_READ_ID_02 = "READ_02";
    public static final String TEST_READ_ID_03 = "READ_03";

    @Test
    public void testSoftclipLeft()
    {
        String refBases = "AAA";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "TTA";
        String cigar1 = "2S1M";

        String readStr2 = "TTA";
        String cigar2 = "2S1M";

        String readStr3 = "GAA";
        String cigar3 = "1S2M";

        String qualStr = DUPLEX_QUAL_STR.repeat(refBases.length());

        SAMRecord read1 = createSbxSamRecord(TEST_READ_ID, CHR_1, 3, readStr1, qualStr, cigar1);
        SAMRecord read2 = createSbxSamRecord(TEST_READ_ID_02, CHR_1, 3, readStr2, qualStr, cigar2);
        SAMRecord read3 = createSbxSamRecord(TEST_READ_ID_03, CHR_1, 2, readStr3, qualStr, cigar3);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "2S1M";
        String expectedReadStr = "TTA";

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(refBases.length(), consensusRead.getAlignmentStart());
        assertEquals(1, consensusRead.getUnclippedStart());
        assertEquals(refBases.length(), consensusRead.getAlignmentEnd());
        assertEquals(refBases.length(), consensusRead.getUnclippedEnd());
    }

    @Test
    public void testSoftclipRight()
    {
        String refBases = "AAA";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "ATT";
        String cigar1 = "1M2S";

        String readStr2 = "ATT";
        String cigar2 = "1M2S";

        String readStr3 = "AAG";
        String cigar3 = "2M1S";

        String qualStr = DUPLEX_QUAL_STR.repeat(refBases.length());

        SAMRecord read1 = createSbxSamRecord(TEST_READ_ID, CHR_1, 1, readStr1, qualStr, cigar1);
        SAMRecord read2 = createSbxSamRecord(TEST_READ_ID_02, CHR_1, 1, readStr2, qualStr, cigar2);
        SAMRecord read3 = createSbxSamRecord(TEST_READ_ID_03, CHR_1, 1, readStr3, qualStr, cigar3);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "1M2S";
        String expectedReadStr = "ATT";

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(1, consensusRead.getAlignmentStart());
        assertEquals(1, consensusRead.getUnclippedStart());
        assertEquals(1, consensusRead.getAlignmentEnd());
        assertEquals(refBases.length(), consensusRead.getUnclippedEnd());
    }

    @Test
    public void testCreateLeftSoftclipOnlyWhenAlignmentScoreImproves()
    {
        String refBases = "ATGC";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "TTGC";
        String cigar = "4M";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord(TEST_READ_ID, CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord(TEST_READ_ID_02, CHR_1, 1, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "4M";

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(1, consensusRead.getAlignmentStart());
        assertEquals(1, consensusRead.getUnclippedStart());
        assertEquals(refBases.length(), consensusRead.getAlignmentEnd());
        assertEquals(refBases.length(), consensusRead.getUnclippedEnd());
    }

    @Test
    public void testCreateRightSoftclipOnlyWhenAlignmentScoreImproves()
    {
        String refBases = "ATGC";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "ATGG";
        String cigar = "4M";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord(TEST_READ_ID, CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord(TEST_READ_ID_02, CHR_1, 1, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "4M";

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(1, consensusRead.getAlignmentStart());
        assertEquals(1, consensusRead.getUnclippedStart());
        assertEquals(refBases.length(), consensusRead.getAlignmentEnd());
        assertEquals(refBases.length(), consensusRead.getUnclippedEnd());
    }
    
    private static SAMRecord createSbxSamRecord(
            final String readName, final String chromosome, int alignmentStart, final String readStr, final String qualStr, final String cigar)
    {
        SAMRecord read = createSamRecordUnpaired(
                readName, chromosome, alignmentStart, readStr, cigar, false, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(60);

        return read;
    }

    private static ConsensusReads getConsensusReads(final String refBases)
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, refBases);
        mockRefGenome.ChromosomeLengths.put(CHR_1, refBases.length());
        return new ConsensusReads(mockRefGenome, SBX);
    }
}
