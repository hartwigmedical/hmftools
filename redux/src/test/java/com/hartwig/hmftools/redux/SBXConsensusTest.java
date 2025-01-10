package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_ERROR_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.ANY_BASE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

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
import com.hartwig.hmftools.redux.consensus.ConsensusOutcome;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.AnnotatedBase;
import com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.ExtendedRefPos;
import com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.SbxBuilder;
import com.hartwig.hmftools.redux.consensus.RefGenome;

import org.junit.Test;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SBXConsensusTest
{
    private final RefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;
    private final SbxBuilder mSbxBuilder;
    private final Map<Byte, Byte> mNextBaseMap;

    public SBXConsensusTest()
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mockRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES.length());

        mConsensusReads = new ConsensusReads(mockRefGenome, SBX);

        mRefGenome = new RefGenome(mockRefGenome);
        mSbxBuilder = new SbxBuilder(mRefGenome);
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

    @Test
    public void testSingleZeroQualBaseInvalidPosition()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A' };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals((byte) 'A', consensus.Base);
        assertEquals((byte) DUPLEX_ERROR_QUAL, consensus.Qual);
        assertEquals(S, consensus.CigarOp);
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testSingleZeroQualBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte[] bases = new byte[] { mNextBaseMap.get(refBase) };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals((byte) 1, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testSingleNoneZeroQualBase()
    {
        int basePosition = 100;
        ExtendedRefPos pos = new ExtendedRefPos(basePosition, 0);
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte[] bases = new byte[] { mNextBaseMap.get(refBase) };
        byte[] quals = new byte[] { (byte) SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(bases[0], consensus.Base);
        assertEquals((byte) SIMPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
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
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals((byte) 1, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
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
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(readBase, consensus.Base);
        assertEquals((byte) SIMPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
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
        byte[] quals = new byte[] { (byte) SIMPLEX_QUAL, (byte) SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals((byte) 1, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
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
        byte[] quals = new byte[] { (byte) DUPLEX_QUAL, (byte) DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(readBase, consensus.Base);
        assertEquals((byte) DUPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
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
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(readBase, consensus.Base);
        assertEquals((byte) SIMPLEX_QUAL, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
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
        byte[] quals = new byte[] { (byte) DUPLEX_QUAL, (byte) DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, M, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(refBase, consensus.Base);
        assertEquals((byte) 1, consensus.Qual);
        assertEquals(M, consensus.CigarOp);
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
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, I, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(ANY_BASE, consensus.Base);
        assertEquals((byte) 0, consensus.Qual);
        assertEquals(I, consensus.CigarOp);
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testZeroQualInvalidPositionSoftClip()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A', (byte) 'A' };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_ERROR_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(ANY_BASE, consensus.Base);
        assertEquals((byte) 0, consensus.Qual);
        assertEquals(S, consensus.CigarOp);
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testNoConsensusSimplexBaseInvalidPositionSoftClip()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A', (byte) 'C' };
        byte[] quals = new byte[] { (byte) SIMPLEX_QUAL, (byte) SIMPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(ANY_BASE, consensus.Base);
        assertEquals((byte) 0, consensus.Qual);
        assertEquals(S, consensus.CigarOp);
        assertTrue(consensus.Annotations.isEmpty());
    }

    @Test
    public void testNoConsensusDuplexBaseInvalidPositionSoftClip()
    {
        ExtendedRefPos pos = new ExtendedRefPos(0, 0);
        byte[] bases = new byte[] { (byte) 'A', (byte) 'C' };
        byte[] quals = new byte[] { (byte) DUPLEX_QUAL, (byte) DUPLEX_QUAL };
        List<AnnotatedBase> annotatedBases = getAnnotatedBases(pos, S, bases, quals);
        AnnotatedBase consensus = mSbxBuilder.determineConsensus(CHR_1, annotatedBases);

        assertEquals(ANY_BASE, consensus.Base);
        assertEquals((byte) 0, consensus.Qual);
        assertEquals(S, consensus.CigarOp);
        assertTrue(consensus.Annotations.isEmpty());
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
}
