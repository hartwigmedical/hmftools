package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_ERROR_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.consensus.BaseBuilderConfig;
import com.hartwig.hmftools.redux.consensus.ConsensusStatistics;
import com.hartwig.hmftools.redux.consensus.RefGenome;
import com.hartwig.hmftools.redux.consensus.SBXBaseBuilderConfig;

import org.junit.Test;

public class SBXConsensusBaseBuildingTest
{
    private final RefGenome mRefGenome;
    private final BaseBuilderConfig mBaseBuilder;
    private final Map<Byte, Byte> mNextBaseMap;

    public SBXConsensusBaseBuildingTest()
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mockRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES.length());
        mRefGenome = new RefGenome(mockRefGenome);
        mBaseBuilder = new SBXBaseBuilderConfig(mRefGenome, new ConsensusStatistics());

        mNextBaseMap = Maps.newHashMap();
        mNextBaseMap.put((byte) 'G', (byte) 'C');
        mNextBaseMap.put((byte) 'C', (byte) 'A');
        mNextBaseMap.put((byte) 'A', (byte) 'T');
        mNextBaseMap.put((byte) 'T', (byte) 'G');
    }

    @Test
    public void testSingleZeroQualBaseInvalidPosition()
    {
        byte[] bases = new byte[] { (byte) 'A' };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, INVALID_POSITION);

        assertEquals(consensusBaseAndQual[0], (byte) 'A');
        assertEquals(consensusBaseAndQual[1], (byte) DUPLEX_ERROR_QUAL);
    }

    @Test
    public void testSingleZeroQualBase()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte[] bases = new byte[] { mNextBaseMap.get(refBase) };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], refBase);
        assertEquals(consensusBaseAndQual[1], (byte) 1);
    }

    @Test
    public void testSingleNoneZeroQualBase()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte[] bases = new byte[] { mNextBaseMap.get(refBase) };
        byte[] quals = new byte[] { (byte) SIMPLEX_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], bases[0]);
        assertEquals(consensusBaseAndQual[1], SIMPLEX_QUAL);
    }

    @Test
    public void testMultipleZeroQualBases()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_ERROR_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], refBase);
        assertEquals(consensusBaseAndQual[1], (byte) 1);
    }

    @Test
    public void testConsensusSimplexBase()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) SIMPLEX_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], readBase);
        assertEquals(consensusBaseAndQual[1], SIMPLEX_QUAL);
    }

    @Test
    public void testNoConsensusSimplexBases()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase1 = mNextBaseMap.get(refBase);
        byte readBase2 = mNextBaseMap.get(readBase1);
        byte[] bases = new byte[] { readBase1, readBase2 };
        byte[] quals = new byte[] { (byte) SIMPLEX_QUAL, (byte) SIMPLEX_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], refBase);
        assertEquals(consensusBaseAndQual[1], (byte) 1);
    }

    @Test
    public void testConsensusDuplexBase()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase };
        byte[] quals = new byte[] { (byte) DUPLEX_QUAL, (byte) DUPLEX_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], readBase);
        assertEquals(consensusBaseAndQual[1], (byte) DUPLEX_QUAL);
    }

    @Test
    public void testDowngradedConsensusDuplexBase()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase = mNextBaseMap.get(refBase);
        byte[] bases = new byte[] { readBase, readBase, readBase };
        byte[] quals = new byte[] { (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_ERROR_QUAL, (byte) DUPLEX_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], readBase);
        assertEquals(consensusBaseAndQual[1], (byte) SIMPLEX_QUAL);
    }

    @Test
    public void testNoConsensusDuplexBase()
    {
        int basePosition = 100;
        byte refBase = mRefGenome.getRefBase(CHR_1, basePosition);
        byte readBase1 = mNextBaseMap.get(refBase);
        byte readBase2 = mNextBaseMap.get(readBase1);
        byte[] bases = new byte[] { readBase1, readBase2 };
        byte[] quals = new byte[] { (byte) DUPLEX_QUAL, (byte) DUPLEX_QUAL };
        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(false, null, bases, quals, CHR_1, basePosition);

        assertEquals(consensusBaseAndQual[0], refBase);
        assertEquals(consensusBaseAndQual[1], (byte) 1);
    }
}
