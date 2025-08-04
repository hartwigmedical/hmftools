package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createConsensusRead;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SbxConsensusTest
{
    private final MockRefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;

    public SbxConsensusTest()
    {
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES.length());
        mConsensusReads = new ConsensusReads(mRefGenome, SequencingType.SBX, new ConsensusStatistics());
    }

    @Test
    public void testSimplexOnlyGroups()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        int position = 10;
        String consensusBases = REF_BASES.substring(position, 20);

        String readBases = consensusBases;

        // reads have a single disagreeing base but find > 50% in agreement
        byte[] baseQuals = buildBaseQuals(readBases.length(), SIMPLEX_QUAL);
        reads.add(createSamRecord(readBases, position, baseQuals));

        readBases = consensusBases.substring(0, 4) + MockRefGenome.getNextBase(consensusBases.substring(4, 5)) + consensusBases.substring(5);
        baseQuals = buildBaseQuals(readBases.length(), SIMPLEX_QUAL);
        reads.add(createSamRecord(readBases, position, baseQuals));

        readBases = consensusBases.substring(0, 6) + MockRefGenome.getNextBase(consensusBases.substring(6, 7)) + consensusBases.substring(7);
        baseQuals = buildBaseQuals(readBases.length(), SIMPLEX_QUAL);
        reads.add(createSamRecord(readBases, position, baseQuals));

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(ALIGNMENT_ONLY, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());

        for(int i = 0; i < readInfo.ConsensusRead.getBaseQualities().length; ++i)
        {
            assertEquals(SIMPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[i]);
        }

        assertEquals(position, readInfo.ConsensusRead.getAlignmentStart());
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
