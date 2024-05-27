package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_PARTITION_SIZE;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_POS_BUFFER_SIZE;

import java.util.List;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.common.Fragment;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public final class TestUtils
{
    public static final String TEST_READ_BASES = MockRefGenome.generateRandomBases(100);
    public static final String TEST_READ_ID = "READ_01";
    public static final String TEST_READ_CIGAR = "100M";

    public static final String REF_BASES_A = "AAAAAAAAAA";
    public static final String REF_BASES_C = "CCCCCCCCCC";
    public static final String REF_BASES_G = "GGGGGGGGGG";
    public static final String REF_BASES_T = "TTTTTTTTTT";

    public static final String REF_BASES_RANDOM = generateRandomBases(10);

    public static final String REF_BASES = "X" + REF_BASES_RANDOM + REF_BASES_A + REF_BASES_C + REF_BASES_G + REF_BASES_T
            + REF_BASES_A + REF_BASES_C + REF_BASES_G + REF_BASES_T + REF_BASES_RANDOM;

    public static final String REF_BASES_REPEAT_40 = REF_BASES.repeat(40);

    public static final int DEFAULT_QUAL = SamRecordTestUtils.DEFAULT_BASE_QUAL;

    public static ReduxConfig createTestConfig()
    {
        return new ReduxConfig(DEFAULT_PARTITION_SIZE, DEFAULT_POS_BUFFER_SIZE, new MockRefGenome(), false, false, false);
    }

    public static Fragment createFragment(final String readId, final String chrStr, int readStart)
    {
        SAMRecord read = createSamRecord(readId, chrStr, readStart, TEST_READ_BASES, TEST_READ_CIGAR, chrStr, 200,
                false, false, null);
        return new Fragment(read);
    }

    public static Fragment createFragment(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar, final String mateChr,
            int mateStart, boolean isReversed, boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        SAMRecord read = createSamRecord(readId, chrStr, readStart, readBases, cigar, mateChr, mateStart,
                isReversed, isSupplementary, suppAlignment);
        return new Fragment(read);
    }

    public static Fragment createFragment(
            final String readId, final String chrStr, int readStart, final String cigar, boolean isReversed,
            final String mateChr, int mateStart, boolean mateReversed, final String mateCigar)
    {
        SAMRecord read = createSamRecord(
                readId, chrStr, readStart, TEST_READ_BASES, cigar, mateChr, mateStart, isReversed, false, null);

        read.setAttribute(MATE_CIGAR_ATTRIBUTE, mateCigar);
        read.setMateNegativeStrandFlag(mateReversed);
        return new Fragment(read);
    }

    public static void setBaseQualities(final Fragment fragment, int value)
    {
        fragment.reads().forEach(x -> setBaseQualities(x, value));
    }

    public static void setBaseQualities(final SAMRecord read, int value)
    {
        for(int i = 0; i < read.getBaseQualities().length; ++i)
            read.getBaseQualities()[i] = (byte)value;
    }

    public static void setSecondInPair(final SAMRecord read)
    {
        read.setFirstOfPairFlag(false);
        read.setSecondOfPairFlag(true);
    }

    public static ConsensusReadInfo createConsensusRead(final ConsensusReads consensusReads, final List<SAMRecord> reads, final String umiId)
    {
        return consensusReads.createConsensusRead(reads, null, null, umiId);
    }
}
