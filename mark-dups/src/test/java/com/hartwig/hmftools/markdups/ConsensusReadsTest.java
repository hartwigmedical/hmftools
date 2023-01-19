package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.ALIGNMENT_ONLY;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.markdups.umi.ConsensusReadInfo;
import com.hartwig.hmftools.markdups.umi.ConsensusReads;
import com.hartwig.hmftools.markdups.umi.UmiConfig;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ConsensusReadsTest
{
    private final MockRefGenome mRefGenome;
    private final UmiConfig mConfig;
    private final ConsensusReads mConsensusReads;
    private final ReadIdGenerator mReadIdGen;

    private final String REF_BASES_A = "AAAAAAAAAA";
    private final String REF_BASES_C = "CCCCCCCCCC";
    private final String REF_BASES_G = "GGGGGGGGGG";
    private final String REF_BASES_T = "TTTTTTTTTT";
    private final String REF_BASES_RANDOM = generateRandomBases(10);

    private final String REF_BASES = "X" + REF_BASES_RANDOM + REF_BASES_A + REF_BASES_C + REF_BASES_G + REF_BASES_T
            + REF_BASES_A + REF_BASES_C + REF_BASES_G + REF_BASES_T + REF_BASES_RANDOM;

    private final String UMI_ID_1 = "TAGTAG";

    public ConsensusReadsTest()
    {
        mConfig = new UmiConfig(true, true);
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mConsensusReads = new ConsensusReads(mConfig, mRefGenome);
        mReadIdGen = new ReadIdGenerator();
    }

    @Test
    public void testBasicConsensusReads()
    {
        final List<SAMRecord> reads = Lists.newArrayList();

        int posStart = 11;
        String consensusBases = REF_BASES.substring(posStart, 21);
        reads.add(createSamRecord(nextReadId(), posStart, consensusBases, "10M", false));
        reads.add(createSamRecord(nextReadId(), posStart, consensusBases, "10M", false));
        reads.add(createSamRecord(nextReadId(), posStart, consensusBases, "10M", false));

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(reads, UMI_ID_1);
        assertEquals(ALIGNMENT_ONLY, readInfo.Outcome);

    }

    private static SAMRecord createSamRecord(
            final String readId, int readStart, final String readBases, final String cigar, boolean isReversed)
    {
        return TestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, 5000, isReversed, false, null);
    }

    private String nextReadId() { return nextReadId(UMI_ID_1); }

    private String nextReadId(final String umiId)
    {
        String readId = mReadIdGen.nextId();
        return format("ABCD:%s:%s", readId, umiId);
    }
}
