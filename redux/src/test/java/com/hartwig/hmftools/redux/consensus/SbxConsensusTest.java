package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.ANY_BASE;
import static com.hartwig.hmftools.redux.consensus.SbxBaseBuilder.DUPLEX_ERROR_QUAL;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.AnnotatedBase;
import com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.ExtendedRefPos;

// import org.junit.jupiter.api.Test;
import org.junit.Ignore;
import org.junit.Test;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SbxConsensusTest
{
    private static final String ZERO_QUAL_STR = String.valueOf(phredToFastq(0));
    private static final String SIMPLEX_QUAL_STR = String.valueOf(phredToFastq(SIMPLEX_QUAL));
    private static final String DUPLEX_QUAL_STR = String.valueOf(phredToFastq(DUPLEX_QUAL));

    private final RefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;
    private final SbxBaseBuilder mSbxBuilder;
    private final Map<Byte, Byte> mNextBaseMap;

    public SbxConsensusTest()
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mockRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES.length());

        mConsensusReads = new ConsensusReads(mockRefGenome, SBX);

        mRefGenome = new RefGenome(mockRefGenome);
        mSbxBuilder = new SbxBaseBuilder(mRefGenome);
        mSbxBuilder.setChromosomeLength(REF_BASES.length());

        mNextBaseMap = Maps.newHashMap();
        mNextBaseMap.put((byte) 'G', (byte) 'C');
        mNextBaseMap.put((byte) 'C', (byte) 'A');
        mNextBaseMap.put((byte) 'A', (byte) 'T');
        mNextBaseMap.put((byte) 'T', (byte) 'G');
    }

    @Test
    public void testSbxConsensusReadSimple()
    {
        int alignmentStart = 20;
        int readLength = 50;
        String readStr = REF_BASES.substring(alignmentStart - 1, alignmentStart - 1 + readLength);
        String qualStr = String.valueOf(phredToFastq(DUPLEX_QUAL)).repeat(readLength);
        String cigar = readLength + "M";

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, alignmentStart, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, alignmentStart, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);

        ConsensusReadInfo consensusOutput = mConsensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(alignmentStart, consensusRead.getAlignmentStart());
        assertEquals(cigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
    }

    @Ignore
    @Test
    public void testSbxConsensusReadWithClippingAndIndels()
    {
        // construct read1 and read2
        String qualStr = String.valueOf(phredToFastq(DUPLEX_QUAL)).repeat(50);
        StringBuilder readStr = new StringBuilder();
        String cigar = "1S10M1D10M1I28M";

        readStr.append((char) ((byte) mNextBaseMap.get(mRefGenome.getRefBase(CHR_1, 20))));
        for(int i = 21; i < 31; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        for(int i = 32; i < 42; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        readStr.append((char) mRefGenome.getRefBase(CHR_1, 42));
        for(int i = 42; i < 70; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        assertEquals(qualStr.length(), readStr.length());
        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 21, readStr.toString(), qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 21, readStr.toString(), qualStr, cigar);

        // construct read3
        readStr = new StringBuilder();
        cigar = "10M1D40M";

        for(int i = 21; i < 31; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        for(int i = 32; i < 72; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        assertEquals(qualStr.length(), readStr.length());
        SAMRecord read3 = createSbxSamRecord("READ_003", CHR_1, 21, readStr.toString(), qualStr, cigar);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = mConsensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 21;
        String expectedCigar = "1S10M1D10M1I30M";
        StringBuilder expectedReadStr = new StringBuilder();
        expectedReadStr.append((char) ((byte) mNextBaseMap.get(mRefGenome.getRefBase(CHR_1, 20))));
        for(int i = 21; i < 31; i++)
            expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        for(int i = 32; i < 42; i++)
            expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, 42));
        for(int i = 42; i < 72; i++)
            expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        String expectedQualStr = String.valueOf(phredToFastq(DUPLEX_QUAL)).repeat(expectedReadStr.length());

        assertEquals(INDEL_MISMATCH, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr.toString(), consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }

    @Ignore
    @Test
    public void testSbxConsensusDropQualZeroInserts()
    {
        String refBases = "AA";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "ATA";
        String readStr2 = "AGA";
        String readStr3 = "ACA";

        String qualStr = DUPLEX_QUAL_STR + SIMPLEX_QUAL_STR + DUPLEX_QUAL_STR;
        String cigar = "1M1I1M";

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr1, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr2, qualStr, cigar);
        SAMRecord read3 = createSbxSamRecord("READ_003", CHR_1, 1, readStr3, qualStr, cigar);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 1;
        String expectedCigar = "2M";
        String expectedReadStr = "AA";
        String expectedQualStr = DUPLEX_QUAL_STR.repeat(2);

        assertEquals(INDEL_MATCH, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }

    @Ignore
    @Test
    public void testSbxConsensusDropNonStrictMajorityInserts()
    {
        String refBases = "A".repeat(20) + "AAAAAAAA" + "A".repeat(20);
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "A".repeat(20) + "AATAATAATAA" + "A".repeat(20);
        String qualStr1 = SIMPLEX_QUAL_STR.repeat(readStr1.length());
        String cigar1 = "22M1I2M1I2M1I22M";

        String readStr2 = "A".repeat(20) + "AAAATAATAA" + "A".repeat(20);
        String qualStr2 = SIMPLEX_QUAL_STR.repeat(readStr2.length());
        String cigar2 = "24M1I2M1I22M";

        String readStr3 = "A".repeat(20) + "AAAAAATAA" + "A".repeat(20);
        String qualStr3 = SIMPLEX_QUAL_STR.repeat(readStr3.length());
        String cigar3 = "26M1I22M";

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr1, qualStr1, cigar1);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr2, qualStr2, cigar2);
        SAMRecord read3 = createSbxSamRecord("READ_003", CHR_1, 1, readStr3, qualStr3, cigar3);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 1;
        String expectedCigar = cigar2;
        String expectedReadStr = readStr2;
        String expectedQualStr = qualStr2;

        assertEquals(INDEL_MISMATCH, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }

    @Ignore
    @Test
    public void testSbxConsensusKeepStrictMajorityDels()
    {
        String refBases = "A".repeat(20) + "ATGCATGC" + "A".repeat(20);
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "A".repeat(20) + "AGCTC" + "A".repeat(20);
        String qualStr1 = SIMPLEX_QUAL_STR.repeat(readStr1.length());
        String cigar1 = "21M1D2M1D1M1D21M";

        String readStr2 = "A".repeat(20) + "ATGCTC" + "A".repeat(20);
        String qualStr2 = SIMPLEX_QUAL_STR.repeat(readStr2.length());
        String cigar2 = "24M1D1M1D21M";

        String readStr3 = "A".repeat(20) + "ATGCATC" + "A".repeat(20);
        String qualStr3 = SIMPLEX_QUAL_STR.repeat(readStr3.length());
        String cigar3 = "26M1D21M";

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr1, qualStr1, cigar1);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr2, qualStr2, cigar2);
        SAMRecord read3 = createSbxSamRecord("READ_003", CHR_1, 1, readStr3, qualStr3, cigar3);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 1;
        String expectedCigar = cigar2;
        String expectedReadStr = readStr2;
        String expectedQualStr = qualStr2;

        assertEquals(INDEL_MISMATCH, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }

    @Ignore
    @Test
    public void testSbxConsensusNoReplacementOfSoftClipWithRef()
    {
        String refBases = "AA";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "TA";
        String readStr2 = "CA";

        String qualStr = SIMPLEX_QUAL_STR.repeat(readStr1.length());
        String cigar = "1S1M";

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 2, readStr1, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 2, readStr2, qualStr, cigar);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 2;
        String expectedCigar = "1S1M";
        String expectedReadStr = "NA";
        String expectedQualStr = ZERO_QUAL_STR + SIMPLEX_QUAL_STR;

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }

    @Test
    public void testSbxConsensusAlignmentBoundariesForInsertBookends()
    {
        int alignmentStart = 50;
        int alignmentLength = 10;

        String readStr = REF_BASES.charAt(alignmentStart - 1) +
                REF_BASES.substring(alignmentStart - 1, alignmentStart - 1 + alignmentLength) +
                REF_BASES.charAt(alignmentStart - 1 + alignmentLength);

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());
        String cigar = "1I" + alignmentLength + "M1I";

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, alignmentStart, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, alignmentStart, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        ConsensusState consensusState = new ConsensusState(true, CHR_1, mRefGenome);
        mSbxBuilder.buildConsensusRead(reads, consensusState, true);

        assertEquals(alignmentStart, consensusState.MinAlignedPosStart);
        assertEquals(alignmentStart + alignmentLength - 1, consensusState.MaxAlignedPosEnd);
        assertEquals(alignmentStart - 1, consensusState.MinUnclippedPosStart);
        assertEquals(alignmentStart + alignmentLength, consensusState.MaxUnclippedPosEnd);
    }

    @Test
    public void testSingleZeroQualBaseInvalidPosition()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A' };
        byte[] quals = new byte[] { DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals((byte) 'A', consensus.Base);
        assertEquals(DUPLEX_ERROR_QUAL, consensus.Qual);
        assertEquals(S, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testSingleZeroQualBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte[] bases = new byte[] { mNextBaseMap.get(refBase) };
        byte[] quals = new byte[] { DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals((byte) 1, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testSingleNoneZeroQualBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte[] bases = new byte[] { mNextBaseMap.get(refBase) };
        byte[] quals = new byte[] { SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(bases[0], consensus.Base);
        assertEquals(SIMPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testMultipleZeroQualBases()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase };
        byte[] quals = new byte[] { DUPLEX_ERROR_QUAL, DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals((byte) 1, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testConsensusSimplexBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase };
        byte[] quals = new byte[] { DUPLEX_ERROR_QUAL, SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(readBase, consensus.Base);
        assertEquals(SIMPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testNoConsensusSimplexBases()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase1 = mNextBaseMap.get(refBase);
        byte readBase2 = mNextBaseMap.get(readBase1);
        byte[] bases = new byte[] { readBase1, readBase2 };
        byte[] quals = new byte[] { SIMPLEX_QUAL, SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals((byte) 1, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testConsensusDuplexBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase };
        byte[] quals = new byte[] { DUPLEX_QUAL, DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(readBase, consensus.Base);
        assertEquals(DUPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testDowngradedConsensusDuplexBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase, readBase };
        byte[] quals = new byte[] { DUPLEX_ERROR_QUAL, DUPLEX_ERROR_QUAL, DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(readBase, consensus.Base);
        assertEquals(SIMPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testNoConsensusDuplexBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase1 = mNextBaseMap.get(refBase);
        byte readBase2 = mNextBaseMap.get(readBase1);
        byte[] bases = new byte[] { readBase1, readBase2 };
        byte[] quals = new byte[] { DUPLEX_QUAL, DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals(1, consensus.Qual);
        assertEquals(M, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testZeroQualInsert()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 1);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase };
        byte[] quals = new byte[] { DUPLEX_ERROR_QUAL, DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, I, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertNull(consensus);
    }

    @Test
    public void testZeroQualInvalidPositionSoftClip()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A', (byte) 'A' };
        byte[] quals = new byte[] { DUPLEX_ERROR_QUAL, DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(ANY_BASE, consensus.Base);
        assertEquals((byte) 0, consensus.Qual);
        assertEquals(S, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testNoConsensusSimplexBaseInvalidPositionSoftClip()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A', (byte) 'C' };
        byte[] quals = new byte[] { SIMPLEX_QUAL, SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(ANY_BASE, consensus.Base);
        assertEquals((byte) 0, consensus.Qual);
        assertEquals(S, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testNoConsensusDuplexBaseInvalidPositionSoftClip()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A', (byte) 'C' };
        byte[] quals = new byte[] { DUPLEX_QUAL, DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(ANY_BASE, consensus.Base);
        assertEquals((byte) 0, consensus.Qual);
        assertEquals(S, consensus.cigarOp());
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testDropMultiMaxDuplexInsert()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 1);
        byte[] bases = new byte[] { (byte) 'A', (byte) 'C' };
        byte[] quals = new byte[] { DUPLEX_QUAL, DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, I, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertNull(consensus);
    }

    @Test
    public void testDropSingleBaseQualZeroInsert()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 1);
        byte[] bases = new byte[] { (byte) 'A' };
        byte[] quals = new byte[] { (byte) 0 };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, I, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertNull(consensus);
    }

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

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 3, readStr1, qualStr, cigar1);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 3, readStr2, qualStr, cigar2);
        SAMRecord read3 = createSbxSamRecord("READ_003", CHR_1, 2, readStr3, qualStr, cigar3);
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

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr1, qualStr, cigar1);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr2, qualStr, cigar2);
        SAMRecord read3 = createSbxSamRecord("READ_003", CHR_1, 1, readStr3, qualStr, cigar3);
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

    @Ignore
    @Test
    public void testInsertToSoftclipLeft()
    {
        String refBases = "A".repeat(3);
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "A".repeat(4);
        String cigar = "1I3M";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "1S3M";

        assertEquals(INDEL_MATCH, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(1, consensusRead.getAlignmentStart());
        assertEquals(0, consensusRead.getUnclippedStart());
        assertEquals(refBases.length(), consensusRead.getAlignmentEnd());
        assertEquals(refBases.length(), consensusRead.getUnclippedEnd());
    }

    @Ignore
    @Test
    public void testInsertToSoftclipRight()
    {
        String refBases = "A".repeat(4);
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "A".repeat(3) + "G";
        String cigar = "3M1I";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "3M1S";

        assertEquals(INDEL_MATCH, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(1, consensusRead.getAlignmentStart());
        assertEquals(1, consensusRead.getUnclippedStart());
        assertEquals(3, consensusRead.getAlignmentEnd());
        assertEquals(readStr.length(), consensusRead.getUnclippedEnd());
    }

    @Ignore
    @Test
    public void testMultiInsertToSoftclipLeft()
    {
        String refBases = "A".repeat(3);
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "A".repeat(5);
        String cigar = "2I3M";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "2S3M";

        assertEquals(INDEL_MATCH, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(1, consensusRead.getAlignmentStart());
        assertEquals(-1, consensusRead.getUnclippedStart());
        assertEquals(refBases.length(), consensusRead.getAlignmentEnd());
        assertEquals(refBases.length(), consensusRead.getUnclippedEnd());
    }

    @Ignore
    @Test
    public void testDeleteAbsorbedIntoSoftclipLeft()
    {
        String refBases = "ATGC";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "AGC";
        String cigar = "1M1D2M";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "1S2M";

        assertEquals(INDEL_MATCH, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(3, consensusRead.getAlignmentStart());
        assertEquals(2, consensusRead.getUnclippedStart());
        assertEquals(refBases.length(), consensusRead.getAlignmentEnd());
        assertEquals(refBases.length(), consensusRead.getUnclippedEnd());
    }

    @Ignore
    @Test
    public void testDeleteAbsorbedIntoSoftclipRight()
    {
        String refBases = "ATGC";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "ATC";
        String cigar = "2M1D1M";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr, qualStr, cigar);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);

        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        String expectedCigar = "2M1S";

        assertEquals(INDEL_MATCH, consensusOutcome);
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(1, consensusRead.getAlignmentStart());
        assertEquals(1, consensusRead.getUnclippedStart());
        assertEquals(2, consensusRead.getAlignmentEnd());
        assertEquals(3, consensusRead.getUnclippedEnd());
    }

    @Test
    public void testCreateLeftSoftclipOnlyWhenAlignmentScoreImproves()
    {
        String refBases = "ATGC";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr = "TTGC";
        String cigar = "4M";

        String qualStr = DUPLEX_QUAL_STR.repeat(readStr.length());

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr, qualStr, cigar);
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

        SAMRecord read1 = createSbxSamRecord("READ_001", CHR_1, 1, readStr, qualStr, cigar);
        SAMRecord read2 = createSbxSamRecord("READ_002", CHR_1, 1, readStr, qualStr, cigar);
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

    private static List<AnnotatedBase> getAnnotatedBases(
            final ExtendedRefPos pos, final CigarOperator cigarOp, final byte[] bases, final byte[] quals)
    {
        List<AnnotatedBase> annotatedBases = Lists.newArrayList();
        for(int i = 0; i < bases.length; i++)
        {
            byte base = bases[i];
            byte qual = quals[i];
            annotatedBases.add(new AnnotatedBase(pos, base, qual, cigarOp));
        }

        return annotatedBases;
    }

    private static SAMRecord createSbxSamRecord(final String readName, final String chromosome, int alignmentStart, final String readStr,
            final String qualStr, final String cigar)
    {
        SAMRecord read = createSamRecordUnpaired(readName, chromosome, alignmentStart, readStr, cigar, false, false, null);
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
