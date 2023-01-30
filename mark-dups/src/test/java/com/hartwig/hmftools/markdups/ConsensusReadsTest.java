package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.markdups.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.markdups.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.markdups.umi.IndelConsensusReads.haveConsistentCigars;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

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
        mConfig = new UmiConfig(true);
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
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());
        assertEquals(posStart, readInfo.ConsensusRead.getAlignmentStart());

        // mismatch at some of the bases and so determined by qual
        reads.clear();

        SAMRecord read1 = createSamRecord(nextReadId(), posStart, consensusBases, "10M", false);
        setBaseQualities(read1, DEFAULT_QUAL - 1);
        reads.add(read1);

        String bases2 = "AATAATAATA";
        SAMRecord read2 = createSamRecord(nextReadId(), posStart + 2, bases2, "2S8M", false);
        reads.add(read2);

        String bases3 = "AATAAGAAAA";
        SAMRecord read3 = createSamRecord(nextReadId(), posStart + 4, bases3, "4S6M", false);
        setBaseQualities(read3, DEFAULT_QUAL - 1);
        reads.add(read3);

        // 3rd base is 'T' due to max qual of 2nd and 3rd reads
        // 6th base is 'T' due to max qual of 2nd read
        // 9th base is 'A' due to max qual of 1st and 3rd reads

        readInfo = mConsensusReads.createConsensusRead(reads, UMI_ID_1);
        assertEquals(ALIGNMENT_ONLY, readInfo.Outcome);
        assertEquals(posStart, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());

        assertEquals("AATAATAAAA", readInfo.ConsensusRead.getReadString());
        assertEquals(19, readInfo.ConsensusRead.getBaseQualities()[2]);
        assertEquals(0, readInfo.ConsensusRead.getBaseQualities()[5]);
        assertEquals(18, readInfo.ConsensusRead.getBaseQualities()[8]);

        // test preference for reference base when quals are qual
        reads.clear();
        posStart = 16;
        read1 = createSamRecord(nextReadId(), posStart, REF_BASES_A, "10M", false);
        reads.add(read1);

        read2 = createSamRecord(nextReadId(), posStart, REF_BASES_C, "10M", false);
        reads.add(read2);

        readInfo = mConsensusReads.createConsensusRead(reads, UMI_ID_1);
        assertEquals("AAAAACCCCC", readInfo.ConsensusRead.getReadString());

        // soft-clips at both ends
        reads.clear();
        posStart = 11;
        read1 = createSamRecord(nextReadId(), posStart + 3, REF_BASES_A, "3S6M1S", false);
        reads.add(read1);

        read2 = createSamRecord(nextReadId(), posStart + 2, REF_BASES_A, "2S5M3S", false);
        reads.add(read2);

        readInfo = mConsensusReads.createConsensusRead(reads, UMI_ID_1);

        assertEquals(posStart + 2, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals("2S7M1S", readInfo.ConsensusRead.getCigarString());
        assertEquals(REF_BASES_A, readInfo.ConsensusRead.getReadString());

        // longer reads at the non-fragment start end
        reads.clear();
        posStart = 11;
        read1 = createSamRecord(nextReadId(), posStart, REF_BASES_A, "8M2S", false);
        reads.add(read1);

        read2 = createSamRecord(nextReadId(), posStart, REF_BASES_A + "CCC", "7M6S", false);
        reads.add(read2);

        readInfo = mConsensusReads.createConsensusRead(reads, UMI_ID_1);

        assertEquals(posStart, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals("8M5S", readInfo.ConsensusRead.getCigarString());
        assertEquals(REF_BASES_A + "CCC", readInfo.ConsensusRead.getReadString());

        // same again on the reverse strand
        reads.clear();
        posStart = 11;
        read1 = createSamRecord(nextReadId(), posStart + 2, REF_BASES_A, "2S6M2S", true);
        reads.add(read1);

        read2 = createSamRecord(nextReadId(), posStart, "CCCC" + REF_BASES_A, "4S7M3S", true);
        reads.add(read2);

        readInfo = mConsensusReads.createConsensusRead(reads, UMI_ID_1);

        assertEquals(posStart, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals("4S8M2S", readInfo.ConsensusRead.getCigarString());
        assertEquals("CCCC" + REF_BASES_A, readInfo.ConsensusRead.getReadString());
    }

    @Test
    public void testIndelCompatibility()
    {
        int posStart = 13;

        String consensusBases = REF_BASES_A.substring(0, 5) + "C" + REF_BASES_A.substring(5, 10);
        String firstCigar = "2S3M1I3M2S";
        SAMRecord read1 = createSamRecord(nextReadId(), posStart, consensusBases, firstCigar, false);

        SAMRecord read2 = createSamRecord(nextReadId(), posStart, consensusBases, firstCigar, false);
        assertTrue(haveConsistentCigars(Lists.newArrayList(read1, read2)));

        // differing soft-clipping
        read2 = createSamRecord(nextReadId(), posStart, consensusBases, "3S3M1I3M3S", false);
        assertTrue(haveConsistentCigars(Lists.newArrayList(read1, read2)));

        // differing initial and end alignment
        read2 = createSamRecord(nextReadId(), posStart, consensusBases, "10M1I5M", false);
        assertTrue(haveConsistentCigars(Lists.newArrayList(read1, read2)));

        // differing initial and end alignment and soft-clipping
        read2 = createSamRecord(nextReadId(), posStart, consensusBases, "3S4M1I4M3S", false);
        assertTrue(haveConsistentCigars(Lists.newArrayList(read1, read2)));

        // now test differences

        // different indel length
        read2 = createSamRecord(nextReadId(), posStart, consensusBases, "2S3M2I3M2S", false);
        assertFalse(haveConsistentCigars(Lists.newArrayList(read1, read2)));

        firstCigar = "5M2D5M1I5M";
        read1 = createSamRecord(nextReadId(), posStart, consensusBases, firstCigar, false);
        read2 = createSamRecord(nextReadId(), posStart, consensusBases, firstCigar, false);
        assertTrue(haveConsistentCigars(Lists.newArrayList(read1, read2)));

        read2 = createSamRecord(nextReadId(), posStart, consensusBases, "5M2D4M1I5M", false);
        assertFalse(haveConsistentCigars(Lists.newArrayList(read1, read2)));
    }

    @Test
    public void testMatchingIndelReads()
    {
        final List<SAMRecord> reads = Lists.newArrayList();

        int posStart = 13;

        String consensusBases = REF_BASES_A.substring(0, 5) + "C" + REF_BASES_A.substring(5, 10);
        String indelCigar = "2S3M1I3M2S";
        SAMRecord read1 = createSamRecord(nextReadId(), posStart, consensusBases, indelCigar, false);
        reads.add(read1);

        SAMRecord read2 = createSamRecord(nextReadId(), posStart, consensusBases, indelCigar, false);
        reads.add(read2);

        SAMRecord read3 = createSamRecord(nextReadId(), posStart, consensusBases, indelCigar, false);
        reads.add(read3);

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(reads, UMI_ID_1);
        assertEquals(INDEL_MATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals(indelCigar, readInfo.ConsensusRead.getCigarString());
        assertEquals(posStart, readInfo.ConsensusRead.getAlignmentStart());
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
