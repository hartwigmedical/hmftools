package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_READ_INDEX_TAG;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.extractDuplexBaseIndex;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.redux.ReduxConstants.SBX_DUPLEX_ADJACENT_1_QUAL;
import static com.hartwig.hmftools.redux.ReduxConstants.SBX_DUPLEX_ADJACENT_2_QUAL;
import static com.hartwig.hmftools.redux.ReduxConstants.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.redux.ReduxConstants.SBX_DUPLEX_QUAL;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createConsensusRead;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.DUPLEX_NO_CONSENSUS_QUAL;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.SIMPLEX_NO_CONSENSUS_QUAL;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SbxConsensusTest
{
    private final MockRefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;

    //                                        12345678901234567890
    private static final String REF_BASES = "XAAAAACCCCCGGGGGTTTTT";

    public SbxConsensusTest()
    {
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES.length());
        mConsensusReads = new ConsensusReads(mRefGenome, SequencingType.SBX, new ConsensusStatistics());
    }

    @Test
    public void testSbxSimplexConsensus()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        int position = 1;
        String refBases = REF_BASES.substring(position, 11);

        // reads have a single disagreeing base but find > 50% in agreement
        byte[] baseQuals = buildBaseQuals(refBases.length(), SIMPLEX_QUAL);
        SAMRecord read1 = createSamRecord(refBases, position, baseQuals);
        reads.add(read1);

        String readBases = "T" + refBases.substring(1);
        baseQuals = buildBaseQuals(readBases.length(), SIMPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases, position, baseQuals);
        reads.add(read2);

        baseQuals = buildBaseQuals(readBases.length(), SIMPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases, position, baseQuals);
        reads.add(read3);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(ALIGNMENT_ONLY, readInfo.Outcome);
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());

        for(int i = 0; i < readInfo.ConsensusRead.getBaseQualities().length; ++i)
        {
            assertEquals(SIMPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[i]);
        }

        assertEquals(position, readInfo.ConsensusRead.getAlignmentStart());

        // test again with no bases agreeing
        readBases = "G" + refBases.substring(1);

        read3 = createSamRecord(readBases, position, baseQuals);
        reads.set(2, read3);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(refBases, readInfo.ConsensusRead.getReadString());
        assertEquals(SIMPLEX_NO_CONSENSUS_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);
    }

    @Test
    public void testSbxDuplexConsensus()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        int position = 1;
        String refBases = REF_BASES.substring(position, 11);

        String readBases = "T" + refBases.substring(1);

        // reads all agree
        byte[] baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);
        SAMRecord read1 = createSamRecord(readBases, position, baseQuals);
        reads.add(read1);

        baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases, position, baseQuals);
        reads.add(read2);

        baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases, position, baseQuals);
        reads.add(read3);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());

        for(int i = 0; i < readInfo.ConsensusRead.getBaseQualities().length; ++i)
        {
            assertEquals(DUPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[i]);
        }

        assertEquals(position, readInfo.ConsensusRead.getAlignmentStart());

        // test again with 2 of 3 agreeing, so a consensus is used
        String readBases3 = "G" + refBases.substring(1);

        read3 = createSamRecord(readBases3, position, baseQuals);
        reads.set(2, read3);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals(DUPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);

        // test that low qual reads are considered for the required percentage but not consensus itself
        baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);
        baseQuals[0] = SBX_DUPLEX_MISMATCH_QUAL;
        SAMRecord read4 = createSamRecord(readBases, position, baseQuals);
        reads.add(read4);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(refBases, readInfo.ConsensusRead.getReadString());
        assertEquals(DUPLEX_NO_CONSENSUS_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);

        // a simplex read is not enough for consensus again
        baseQuals = buildBaseQuals(readBases.length(), SIMPLEX_QUAL);
        SAMRecord read5 = createSamRecord(readBases, position, baseQuals);
        reads.add(read5);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(refBases, readInfo.ConsensusRead.getReadString());
        assertEquals(DUPLEX_NO_CONSENSUS_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);

        baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);
        SAMRecord read6 = createSamRecord(readBases, position, baseQuals);
        reads.add(read6);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals(DUPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);
    }

    @Test
    public void testSbxFinaliseReads()
    {
        int position = 1;
        String refBases = REF_BASES.substring(position, 11);

        String readBases = "T" + refBases.substring(1);

        byte[] baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);
        SAMRecord read = createSamRecord(readBases, position, baseQuals);

        SbxRoutines.finaliseRead(mRefGenome, read);

        for(int i = 0; i < read.getBaseQualities().length; ++i)
        {
            assertEquals(SBX_DUPLEX_QUAL, read.getBaseQualities()[i]);
        }

        ConsensusType consensusType = extractConsensusType(read);
        int duplexBaseIndex = SbxRoutines.findMaxDuplexBaseIndex(List.of(read));
        assertEquals(ConsensusType.NONE, consensusType);
        assertEquals(0, duplexBaseIndex);

        baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);
        read = createSamRecord(readBases, position, baseQuals);
        read.setReadNegativeStrandFlag(true);
        read.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, ConsensusType.SINGLE.toString());
        setDupluxBaseIndex(read);

        SbxRoutines.finaliseRead(mRefGenome, read);

        consensusType = extractConsensusType(read);
        duplexBaseIndex = extractDuplexBaseIndex(read);
        assertEquals(ConsensusType.DUAL, consensusType);
        assertEquals(9, duplexBaseIndex);

        // mark the transition point mid-way through the read
        baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);

        int lastIndex = readBases.length() - 1;
        baseQuals[lastIndex] = SIMPLEX_QUAL;
        baseQuals[lastIndex - 1] = SIMPLEX_QUAL;
        baseQuals[lastIndex - 2] = SIMPLEX_QUAL;

        read = createSamRecord(readBases, position, baseQuals);
        read.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, ConsensusType.SINGLE.toString());
        read.setReadNegativeStrandFlag(true);
        setDupluxBaseIndex(read);

        SbxRoutines.finaliseRead(mRefGenome, read);

        consensusType = extractConsensusType(read);
        duplexBaseIndex = extractDuplexBaseIndex(read);
        assertEquals(ConsensusType.DUAL, consensusType);
        assertEquals(6, duplexBaseIndex);

        // test duplex mismatch bases
        baseQuals = buildBaseQuals(readBases.length(), DUPLEX_QUAL);

        baseQuals[3] = SBX_DUPLEX_MISMATCH_QUAL;
        baseQuals[7] = SBX_DUPLEX_MISMATCH_QUAL;
        baseQuals[8] = SBX_DUPLEX_MISMATCH_QUAL;

        read = createSamRecord(readBases, position, baseQuals);
        SbxRoutines.finaliseRead(mRefGenome, read);

        assertEquals(SBX_DUPLEX_QUAL, read.getBaseQualities()[0]);
        assertEquals(SBX_DUPLEX_ADJACENT_2_QUAL, read.getBaseQualities()[1]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL, read.getBaseQualities()[2]);
        assertEquals(SBX_DUPLEX_MISMATCH_QUAL, read.getBaseQualities()[3]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL, read.getBaseQualities()[4]);
        assertEquals(SBX_DUPLEX_ADJACENT_2_QUAL, read.getBaseQualities()[5]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL, read.getBaseQualities()[6]);
        assertEquals(SBX_DUPLEX_MISMATCH_QUAL, read.getBaseQualities()[7]);
        assertEquals(SBX_DUPLEX_MISMATCH_QUAL, read.getBaseQualities()[8]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL, read.getBaseQualities()[9]);
    }

    private static void setDupluxBaseIndex(final SAMRecord record)
    {
        int firstDuplexBaseIndex = SbxRoutines.findMaxDuplexBaseIndex(List.of(record));

        if(firstDuplexBaseIndex >= 0)
            record.setAttribute(SBX_DUPLEX_READ_INDEX_TAG, firstDuplexBaseIndex);
    }

    private static SAMRecord createSamRecord(final String readBases, final int position, final byte[] baseQualities)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, position, readBases, format("%dM", readBases.length()), CHR_1, 300,
                false, false, null, true, TEST_READ_CIGAR);
        record.setBaseQualities(baseQualities);
        return record;
    }
}
