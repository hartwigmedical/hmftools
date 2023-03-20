package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.markdups.ConsensusReadsTest.UMI_ID_1;
import static com.hartwig.hmftools.markdups.ConsensusReadsTest.nextUmiReadId;
import static com.hartwig.hmftools.markdups.TestUtils.REF_BASES;
import static com.hartwig.hmftools.markdups.TestUtils.REF_BASES_A;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.markdups.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.markdups.umi.IndelConsensusReads.haveConsistentCigars;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.markdups.umi.ConsensusReadInfo;
import com.hartwig.hmftools.markdups.umi.ConsensusReads;
import com.hartwig.hmftools.markdups.umi.UmiConfig;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class IndelConsensusReadsTest
{
    private final MockRefGenome mRefGenome;
    private final UmiConfig mConfig;
    private final ConsensusReads mConsensusReads;
    private final ReadIdGenerator mReadIdGen;

    public IndelConsensusReadsTest()
    {
        mConfig = new UmiConfig(true, true);
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mConsensusReads = new ConsensusReads(mConfig, mRefGenome);
        mReadIdGen = new ReadIdGenerator();
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

        read1 = createSamRecord(nextReadId(), posStart, consensusBases, "10M", false);
        read2 = createSamRecord(nextReadId(), posStart, consensusBases, "4M1I6M", false);
        assertFalse(haveConsistentCigars(Lists.newArrayList(read1, read2)));
    }

    @Test
    public void testInsertMismatches()
    {
        // first 2 reads without an insert vs 1 with one ignored
        boolean readsReversed = false;
        SAMRecord read1 = createSamRecord(nextReadId(), 13, REF_BASES_A.substring(0, 10), "2S8M", readsReversed);

        // matching but without the soft-clip
        SAMRecord read2 = createSamRecord(nextReadId(), 11, REF_BASES_A.substring(0, 10), "7M3S", readsReversed);

        String indelBases = REF_BASES_A.substring(0, 5) + "C" + REF_BASES_A.substring(5, 10);
        SAMRecord read3 = createSamRecord(nextReadId(), 12, indelBases, "1S4M1I4M1S", readsReversed);

        indelBases = REF_BASES_A.substring(0, 3) + "GG" + REF_BASES_A.substring(3, 10);
        SAMRecord read4 = createSamRecord(nextReadId(), 11, indelBases, "3M2I7M", readsReversed);

        indelBases = REF_BASES_A.substring(0, 4) + "TTT" + REF_BASES_A.substring(4, 7) + "TTTT" + REF_BASES_A.substring(7, 10);
        SAMRecord read5 = createSamRecord(nextReadId(), 12, indelBases, "1S3M3I3M4I2M1S", readsReversed);

        String consensusBases = REF_BASES_A.substring(0, 10);

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(Lists.newArrayList(read1, read2, read3, read4, read5), UMI_ID_1);
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());
        assertEquals(11, readInfo.ConsensusRead.getAlignmentStart());

        // now with 2 inserts in the consensus
        consensusBases = REF_BASES_A.substring(0, 3) + "TT" + REF_BASES_A.substring(3, 6) + "GGG" + REF_BASES_A.substring(6, 10);
        String consensusCigar = "1S2M2I3M3I3M1S";
        readsReversed = true;
        read1 = createSamRecord(nextReadId(), 12, consensusBases, consensusCigar, readsReversed);

        // matching but without the soft-clip
        read2 = createSamRecord(nextReadId(), 12, consensusBases, consensusCigar, readsReversed);

        read3 = createSamRecord(nextReadId(), 12, REF_BASES_A.substring(0, 10), "2S6M2S", readsReversed);

        // has the same indels but shorter
        indelBases = REF_BASES_A.substring(0, 3) + "TTT" + REF_BASES_A.substring(3, 6) + "GG" + REF_BASES_A.substring(6, 10);
        read4 = createSamRecord(nextReadId(), 12, indelBases, "1S2M3I3M2I3M1S", readsReversed);

        readInfo = mConsensusReads.createConsensusRead(Lists.newArrayList(read1, read2, read3, read4), UMI_ID_1);
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals(consensusCigar, readInfo.ConsensusRead.getCigarString());
        assertEquals(12, readInfo.ConsensusRead.getAlignmentStart());
    }

    @Test
    public void testDeleteMismatches()
    {
        // first 2 reads without a delete vs 1 with one ignored
        SAMRecord read1 = createSamRecord(nextReadId(), 13, REF_BASES_A.substring(0, 10), "2S8M", false);

        // matching but without the soft-clip
        SAMRecord read2 = createSamRecord(nextReadId(), 11, REF_BASES_A.substring(0, 10), "7M3S", false);

        String indelBases = REF_BASES_A.substring(0, 5) + REF_BASES_A.substring(6, 10);
        SAMRecord read3 = createSamRecord(nextReadId(), 12, indelBases, "1S4M1D3M1S", false);

        indelBases = REF_BASES_A.substring(0, 3) + REF_BASES_A.substring(7, 10);
        SAMRecord read4 = createSamRecord(nextReadId(), 11, indelBases, "3M4D3M", false);

        indelBases = REF_BASES_A.substring(0, 3) + REF_BASES_A.substring(4, 6) + REF_BASES_A.substring(7, 10);
        SAMRecord read5 = createSamRecord(nextReadId(), 11, indelBases, "3M1D2M1D3M", false);

        String consensusBases = REF_BASES_A.substring(0, 10);

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(Lists.newArrayList(read1, read2, read3, read4, read5), UMI_ID_1);
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());
        assertEquals(11, readInfo.ConsensusRead.getAlignmentStart());

        // now test a delete in the consensus read
        String consensusCigar = "4M3D3M";
        String delBases = REF_BASES_A.substring(0, 4) + REF_BASES_A.substring(7, 10);
        read1 = createSamRecord(nextReadId(), 11, delBases, consensusCigar, false);
        setBaseQualities(read1, 11);

        // matching but without the soft-clip
        read2 = createSamRecord(nextReadId(), 11, delBases, consensusCigar, false);
        setBaseQualities(read2, 11);

        indelBases = REF_BASES_A.substring(0, 7) + "CCC"; // no delete and differing aligned bases after the delete
        read3 = createSamRecord(nextReadId(), 11, indelBases, "10M", false);
        setBaseQualities(read3, 35);

        consensusBases = REF_BASES_A.substring(0, 4) + "CCC";

        readInfo = mConsensusReads.createConsensusRead(Lists.newArrayList(read1, read2, read3), UMI_ID_1);
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals(consensusCigar, readInfo.ConsensusRead.getCigarString());
        int calcBaseQual = 13;
        assertEquals(calcBaseQual, readInfo.ConsensusRead.getBaseQualities()[4]);
        assertEquals(calcBaseQual, readInfo.ConsensusRead.getBaseQualities()[5]);
        assertEquals(calcBaseQual, readInfo.ConsensusRead.getBaseQualities()[6]);
    }

    @Test
    public void testMismatchedReads()
    {
        // one of the As at the end got deleted
        SAMRecord read1 = createSamRecord(nextReadId(), 10, "CTTCGATAAT", "7M1D3M", true);
        SAMRecord read2 = createSamRecord(nextReadId(), 11, "TTCGATAAAT", "10M", true);

        // they are marked as same read cause read is reversed and the end are the same
        assertEquals(read1.getAlignmentEnd(), read2.getAlignmentEnd());

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(List.of(read1, read2), UMI_ID_1);

        assertEquals("11M", readInfo.ConsensusRead.getCigarString());
        assertEquals("CTTCGATAAAT", readInfo.ConsensusRead.getReadString());
    }

    @Test
    public void testMismatchedReads2()
    {
        // one of the As at the end got deleted
        SAMRecord read1 = createSamRecord(nextReadId(), 10, "CTTCGATAAT", "7M1D3M", true);
        SAMRecord read2 = createSamRecord(nextReadId(), 11, "TTCGATAAAT", "10M", true);

        // they are marked as same read cause read is reversed and the end are the same
        assertEquals(read1.getAlignmentEnd(), read2.getAlignmentEnd());

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(List.of(read1, read2), UMI_ID_1);

        assertEquals("11M", readInfo.ConsensusRead.getCigarString());
        assertEquals("CTTCGATAAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());

        SAMRecord read3 = createSamRecord(nextReadId(), 11, "TTCGATAAAT", "10M", true);
        readInfo = mConsensusReads.createConsensusRead(List.of(read1, read2, read3), UMI_ID_1);

        assertEquals("11M", readInfo.ConsensusRead.getCigarString());
        assertEquals("CTTCGATAAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());
    }

    @Test
    public void testMismatchedReadsInsert()
    {
        // one of the As at the end got deleted
        SAMRecord read1 = createSamRecord(nextReadId(), 10, "CTTCGATAAAAT", "7M1I4M", true);
        SAMRecord read2 = createSamRecord(nextReadId(), 11, "TTCGATAAAT", "10M", true);

        // they are marked as same read cause read is reversed and the end are the same
        assertEquals(read1.getAlignmentEnd(), read2.getAlignmentEnd());

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(List.of(read1, read2), UMI_ID_1);

        assertEquals("11M", readInfo.ConsensusRead.getCigarString());
        assertEquals("CTTCGATAAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());

        // add one more read so the first C will get chopped
        SAMRecord read3 = createSamRecord(nextReadId(), 11, "TTCGATAAAT", "10M", true);

        readInfo = mConsensusReads.createConsensusRead(List.of(read1, read2, read3), UMI_ID_1);

        assertEquals("11M", readInfo.ConsensusRead.getCigarString());
        assertEquals("CTTCGATAAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());
    }

    @Test
    public void testSoftClipDeletes()
    {
        // one of the As at the end got deleted
        SAMRecord read1 = createSamRecord(nextReadId(), 40, REF_BASES.substring(40, 100), "30S30M", true);

        SAMRecord read2 = createSamRecord(
                nextReadId(), 7,
                REF_BASES.substring(37, 67) + REF_BASES.substring(70, 100), "30M3D30M", true);

        SAMRecord read3 = createSamRecord(nextReadId(), 40, REF_BASES.substring(40, 100), "30S30M", true);
        SAMRecord read4 = createSamRecord(nextReadId(), 39, REF_BASES.substring(40, 100), "29S31M", true);

        ConsensusReadInfo readInfo = mConsensusReads.createConsensusRead(List.of(read1, read2, read3, read4), UMI_ID_1);

        assertEquals("30M2D31M", readInfo.ConsensusRead.getCigarString());
        String consensusBases = REF_BASES.substring(37, 67) + REF_BASES.substring(69, 100);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals(7, readInfo.ConsensusRead.getAlignmentStart());


        read1 = createSamRecord(
                nextReadId(), 8,
                REF_BASES.substring(8, 38) + REF_BASES.substring(40, 70), "30M2D30M", true);

        read2 = createSamRecord(nextReadId(), 40, REF_BASES.substring(10, 70), "30S30M", true);
        read3 = createSamRecord(nextReadId(), 39, REF_BASES.substring(10, 70), "29S31M", true);

        readInfo = mConsensusReads.createConsensusRead(List.of(read1, read2, read3), UMI_ID_1);

        assertEquals("30M1D31M", readInfo.ConsensusRead.getCigarString());
        consensusBases = REF_BASES.substring(8, 38) + REF_BASES.substring(39, 70);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals(8, readInfo.ConsensusRead.getAlignmentStart());



        /*
        mAlignmentStart = 38446161
        mAlignmentEnd = 38446313
        mMappingQuality = 60
        mCigarString = null
        mCigar = {Cigar@3199} "66M2D85M"

        mAlignmentStart = 38446229
        mAlignmentEnd = 38446313
        mMappingQuality = 60
        mCigarString = null
        mCigar = {Cigar@3211} "66S85M"

        mAlignmentStart = 38446228
        mAlignmentEnd = 38446313
        mMappingQuality = 60
        mCigarString = null
        mCigar = {Cigar@3220} "65S86M"
         */



    }

    private String nextReadId() { return nextUmiReadId(UMI_ID_1, mReadIdGen); }

    private static SAMRecord createSamRecord(
            final String readId, int readStart, final String readBases, final String cigar, boolean isReversed)
    {
        return SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, 5000, isReversed, false, null);
    }

}
