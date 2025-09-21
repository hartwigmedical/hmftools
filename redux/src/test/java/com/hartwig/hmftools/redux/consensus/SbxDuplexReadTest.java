package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_YC_TAG;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndelIndices;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.reverseDuplexIndelIndices;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.cloneSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.redux.ReduxConstants.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.getAnnotatedBases;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.processAnnotatedBases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.common.FragmentCoords;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.CigarOperator;
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
    public void testGetAnnotatedBasesSimple()
    {
        int readLength = 100;
        int alignmentStart = 25;
        String readStr = "A".repeat(readLength);
        String cigar = readLength + "M";
        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        List<Integer> duplexIndelIndices = Collections.emptyList();
        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndelIndices);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);

            assertEquals(i, annotatedBase.ReadIndex);
            assertEquals(i + alignmentStart, annotatedBase.RefPos);
            assertEquals(M, annotatedBase.Op);
            assertEquals('A', (char) annotatedBase.ReadBase);
            assertFalse(annotatedBase.IsDuplexIndel);
        }
    }

    @Test
    public void testGetAnnotatedBasesLeftSoftClip()
    {
        int readLength = 100;
        int alignmentStart = 25;
        int leftSoftClipLength = 10;
        String readStr = "A".repeat(readLength);
        String cigar = format("%dS%dM", leftSoftClipLength, readLength - leftSoftClipLength);
        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        List<Integer> duplexIndelIndices = Collections.emptyList();
        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndelIndices);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);

            CigarOperator expectedOp = i < leftSoftClipLength ? S : M;

            assertEquals(i, annotatedBase.ReadIndex);
            assertEquals(i + alignmentStart - leftSoftClipLength, annotatedBase.RefPos);
            assertEquals(expectedOp, annotatedBase.Op);
            assertEquals('A', (char) annotatedBase.ReadBase);
            assertFalse(annotatedBase.IsDuplexIndel);
        }
    }

    @Test
    public void testGetAnnotatedBasesIndels()
    {
        int readLength = 100;
        int alignmentStart = 25;
        String readStr = "A".repeat(readLength);
        String cigar = "9M1I40M1D50M";
        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        List<Integer> duplexIndelIndices = Collections.emptyList();
        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndelIndices);

        assertEquals(readLength + 1, annotatedBases.size());

        int expectedReadIndex = 0;
        int expectedRefPos = alignmentStart;
        for(int i = 0; i < annotatedBases.size(); i++)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);

            byte expectedReadBase = i != 50 ? (byte) 'A' : INVALID_BASE_QUAL;
            CigarOperator expectedOp = M;
            if(i == 9)
            {
                expectedOp = I;
            }
            else if(i == 50)
            {
                expectedOp = D;
            }

            assertEquals(expectedReadIndex, annotatedBase.ReadIndex);
            assertEquals(expectedRefPos, annotatedBase.RefPos);
            assertEquals(expectedOp, annotatedBase.Op);
            assertEquals(expectedReadBase, annotatedBase.ReadBase);
            assertFalse(annotatedBase.IsDuplexIndel);

            if(i != 49)
            {
                expectedReadIndex++;
            }

            if(i != 8)
            {
                expectedRefPos++;
            }
        }
    }

    @Test
    public void testProcessAnnotatedBasesForwardRead()
    {
        int alignmentStart = 25;

        String readStr = "A".repeat(5) + "CT".repeat(3) + "A".repeat(5);
        String cigar = "5M4I7M";  // inserts are left aligned.
        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);

        String qualStr = NON_ZERO_QUAL.repeat(9) + MISMATCH_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);
        read.setBaseQualityString(qualStr);

        List<Integer> duplexIndelIndices = getDuplexIndelIndices("0-9ZZ5-0");

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(100);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        List<SbxAnnotatedBase> annotedBases = getAnnotatedBases(read, duplexIndelIndices);
        boolean readModified = processAnnotatedBases(refGenome, CHR_1, annotedBases, true);

        List<SbxAnnotatedBase> expectedBases = getAnnotatedBases(read, duplexIndelIndices);

        expectedBases.get(5).setQual(SBX_DUPLEX_MISMATCH_QUAL);
        expectedBases.get(6).setQual(SBX_DUPLEX_MISMATCH_QUAL);

        expectedBases.get(5).deleteBase();
        expectedBases.get(6).deleteBase();

        expectedBases.get(7).setQual(SBX_DUPLEX_MISMATCH_QUAL);
        expectedBases.get(8).setQual(SBX_DUPLEX_MISMATCH_QUAL);

        assertTrue(readModified);
        assertEquals(expectedBases, annotedBases);

        SbxRoutines.SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL = true;
    }

    @Test
    public void testProcessAnnotatedBasesReverseRead()
    {
        int alignmentStart = 25;

        String readStr = "A".repeat(5) + "CT".repeat(3) + "A".repeat(5);
        String cigar = "5M4I7M";  // inserts are left aligned.
        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, true, false, null);

        String qualStr = NON_ZERO_QUAL.repeat(5) + MISMATCH_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(9);
        read.setBaseQualityString(qualStr);

        List<Integer> duplexIndelIndices = getDuplexIndelIndices("0-9ZZ5-0");

        reverseDuplexIndelIndices(duplexIndelIndices, read.getReadBases().length);

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(100);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        List<SbxAnnotatedBase> annotedBases = getAnnotatedBases(read, duplexIndelIndices);
        boolean readModified = processAnnotatedBases(refGenome, CHR_1, annotedBases, false);

        List<SbxAnnotatedBase> expectedBases = getAnnotatedBases(read, duplexIndelIndices);

        expectedBases.get(9).setQual(SBX_DUPLEX_MISMATCH_QUAL);
        expectedBases.get(10).setQual(SBX_DUPLEX_MISMATCH_QUAL);

        expectedBases.get(5).deleteBase();
        expectedBases.get(6).deleteBase();

        expectedBases.get(7).setQual(SBX_DUPLEX_MISMATCH_QUAL);
        expectedBases.get(8).setQual(SBX_DUPLEX_MISMATCH_QUAL);

        assertTrue(readModified);
        assertEquals(expectedBases, annotedBases);
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
        String cigar = "5M2I7M";  // inserts are right-aligned
        // String cigar = "5M2I7M";
        String qualStr = NON_ZERO_QUAL.repeat(7) + MISMATCH_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);
        String ycTagStr = "0-7ZZ5-0";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readStr, cigar, false, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(mapq);
        read.setAttribute(SBX_YC_TAG, ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, nm);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignmentScore);

        SbxRoutines.stripDuplexIndels(refGenome, read);

        SAMRecord newRead = cloneSamRecord(read, READ_ID_GEN.nextId());
        SbxRoutines.stripDuplexIndelsNew(refGenome, newRead);

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

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(100);
        RefGenomeInterface refGenome = getRefGenome(refBases);

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

        SbxRoutines.stripDuplexIndels(refGenome, read);

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

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(100);
        // RefGenomeInterface refGenome = getRefGenome(refBases);

        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, refBases);
        mockRefGenome.RefGenomeMap.put(CHR_2, refBases);
        mockRefGenome.ChromosomeLengths.put(CHR_1, refBases.length());
        mockRefGenome.ChromosomeLengths.put(CHR_2, refBases.length());
        RefGenomeInterface refGenome = new CachedRefGenome(mockRefGenome);

        String alignedBases = REF_BASES.substring(1, 21);
        String suppBases = REF_BASES.substring(30, 54); // "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        String cigar = "20M24S";
        String suppCigar = "20S10M2I12M";
        int suppPosStart = 25;
        String readBases = alignedBases + suppBases;

        SupplementaryReadData suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_POS_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        String ycTagStr = "20-10Z13-0";

        SAMRecord read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);
        read.setAttribute(SBX_YC_TAG, ycTagStr);

        SbxRoutines.stripDuplexIndels(refGenome, read);

        assertEquals("20M23S", read.getCigarString());
        assertEquals(43, read.getReadBases().length);

        // with supplementary having a different aligned and soft-clip length, requiring adjustment
        suppCigar = "18S12M2I12M";

        suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_POS_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);
        read.setAttribute(SBX_YC_TAG, ycTagStr);

        SbxRoutines.stripDuplexIndels(refGenome, read);

        assertEquals("20M23S", read.getCigarString());
        assertEquals(43, read.getReadBases().length);

        // reversed supp data
        suppCigar = "12M2I08M22S";
        suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_NEG_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);
        read.setAttribute(SBX_YC_TAG, ycTagStr);

        SbxRoutines.stripDuplexIndels(refGenome, read);

        assertEquals("20M23S", read.getCigarString());
        assertEquals(43, read.getReadBases().length);

        // repeat with soft-clip on the other side
        cigar = "24S20M";
        suppCigar = "10M2I12M20S";
        readBases = suppBases + alignedBases;

        suppData = new SupplementaryReadData(CHR_2, suppPosStart, SUPP_POS_STRAND, suppCigar, DEFAULT_MAP_QUAL, 0);

        ycTagStr = "10-5Z18-0";

        read = createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, alignmentStart, readBases, cigar, false, false, suppData);
        read.setAttribute(SBX_YC_TAG, ycTagStr);

        SbxRoutines.stripDuplexIndels(refGenome, read);

        assertEquals("23S20M", read.getCigarString());
        assertEquals(43, read.getReadBases().length);
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
