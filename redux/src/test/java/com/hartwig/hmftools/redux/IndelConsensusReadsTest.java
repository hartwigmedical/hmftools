package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.ConsensusReadsTest.UMI_ID_1;
import static com.hartwig.hmftools.redux.ConsensusReadsTest.nextUmiReadId;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_A;
import static com.hartwig.hmftools.redux.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.redux.TestUtils.createConsensusRead;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class IndelConsensusReadsTest
{
    private final MockRefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;
    private final ReadIdGenerator mReadIdGen;

    public IndelConsensusReadsTest()
    {
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mConsensusReads = new ConsensusReads(mRefGenome);
        mReadIdGen = new ReadIdGenerator();
    }

    @Test
    public void testInsertMismatches()
    {
        // first 2 reads without an insert vs 1 with one ignored
        boolean readsReversed = false;
        SAMRecord read1 = createSamRecord(nextReadId(), 11, REF_BASES_A.substring(0, 10), "10M", readsReversed);

        // matching but without the soft-clip
        SAMRecord read2 = createSamRecord(nextReadId(), 11, REF_BASES_A.substring(0, 10), "10M", readsReversed);

        String indelBases = REF_BASES_A.substring(0, 5) + "C" + REF_BASES_A.substring(5, 10);
        SAMRecord read3 = createSamRecord(nextReadId(), 12, indelBases, "1S4M1I4M1S", readsReversed);

        indelBases = REF_BASES_A.substring(0, 3) + "GG" + REF_BASES_A.substring(3, 10);
        SAMRecord read4 = createSamRecord(nextReadId(), 11, indelBases, "3M2I7M", readsReversed);

        indelBases = REF_BASES_A.substring(0, 4) + "TTT" + REF_BASES_A.substring(4, 7) + "TTTT" + REF_BASES_A.substring(7, 10);
        SAMRecord read5 = createSamRecord(nextReadId(), 12, indelBases, "1S3M3I3M4I2M1S", readsReversed);

        String consensusBases = REF_BASES_A.substring(0, 10);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, Lists.newArrayList(read1, read2, read3, read4, read5), UMI_ID_1);
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

        readInfo = createConsensusRead(mConsensusReads, Lists.newArrayList(read1, read2, read3, read4), UMI_ID_1);
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals(consensusCigar, readInfo.ConsensusRead.getCigarString());
        assertEquals(12, readInfo.ConsensusRead.getAlignmentStart());
    }

    @Test
    public void testDeleteMismatches()
    {
        // first 2 reads without a delete vs 1 with one ignored
        SAMRecord read1 = createSamRecord(nextReadId(), 11, REF_BASES_A.substring(0, 10), "10M", false);
        SAMRecord read2 = createSamRecord(nextReadId(), 11, REF_BASES_A.substring(0, 10), "10M", false);

        String indelBases = REF_BASES_A.substring(0, 5) + REF_BASES_A.substring(6, 10);
        SAMRecord read3 = createSamRecord(nextReadId(), 12, indelBases, "1S4M1D3M1S", false);

        indelBases = REF_BASES_A.substring(0, 3) + REF_BASES_A.substring(7, 10);
        SAMRecord read4 = createSamRecord(nextReadId(), 11, indelBases, "3M4D3M", false);

        indelBases = REF_BASES_A.substring(0, 3) + REF_BASES_A.substring(4, 6) + REF_BASES_A.substring(7, 10);
        SAMRecord read5 = createSamRecord(nextReadId(), 11, indelBases, "3M1D2M1D3M", false);

        String consensusBases = REF_BASES_A.substring(0, 10);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, Lists.newArrayList(read1, read2, read3, read4, read5), UMI_ID_1);
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());
        assertEquals(11, readInfo.ConsensusRead.getAlignmentStart());

        // now test a delete in the consensus read
        String consensusCigar = "4M3D3M";
        String delBases = REF_BASES_A.substring(0, 4) + REF_BASES_A.substring(7, 10);
        read1 = createSamRecord(nextReadId(), 11, delBases, consensusCigar, false);
        setBaseQualities(read1, 11);

        read2 = createSamRecord(nextReadId(), 11, delBases, consensusCigar, false);
        setBaseQualities(read2, 11);

        indelBases = REF_BASES_A.substring(0, 7) + "CCC"; // no delete and differing aligned bases after the delete
        read3 = createSamRecord(nextReadId(), 11, indelBases, "10M", false);
        setBaseQualities(read3, 35);

        consensusBases = REF_BASES_A.substring(0, 4) + "CCC";

        readInfo = createConsensusRead(mConsensusReads, Lists.newArrayList(read1, read2, read3), UMI_ID_1);
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
        assertEquals(consensusCigar, readInfo.ConsensusRead.getCigarString());
        int calcBaseQual = 11; // lowered from 13
        assertEquals(calcBaseQual, readInfo.ConsensusRead.getBaseQualities()[4]);
        assertEquals(calcBaseQual, readInfo.ConsensusRead.getBaseQualities()[5]);
        assertEquals(calcBaseQual, readInfo.ConsensusRead.getBaseQualities()[6]);

        // a basic test for following bases correctly
        String readBases = "CTTCGATAATGGCCGGGCCG";
        String shortedReadBases = readBases.substring(0, readBases.length() - 2);
        String shorterClipped = "3S7M10D5M3S";
        read1 = createSamRecord(nextReadId(), 10, shortedReadBases, shorterClipped, false);
        read2 = createSamRecord(nextReadId(), 10, readBases, "5S5M10D5M5S", false);

        readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2), UMI_ID_1);

        assertEquals(shorterClipped, readInfo.ConsensusRead.getCigarString());
        assertEquals(shortedReadBases, readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());
    }

    @Test
    public void testMismatchedReads2()
    {
        // favour the read with the least number of soft-clip bases
        SAMRecord read1 = createSamRecord(nextReadId(), 10, "CTTCGATAAT", "7M1D3M", true);
        SAMRecord read2 = createSamRecord(nextReadId(), 13, "TTCGATAAAT", "2S8M", true);

        // they are marked as same read cause read is reversed and the end are the same
        assertEquals(read1.getAlignmentEnd(), read2.getAlignmentEnd());

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2), UMI_ID_1);

        assertEquals("7M1D3M", readInfo.ConsensusRead.getCigarString());
        assertEquals("CTTCGATAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());

        // a second supporting read results in the extra base being dropped
        SAMRecord read3 = createSamRecord(nextReadId(), 13, "TTCGATAAAT", "2S8M", true);
        readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2, read3), UMI_ID_1);

        assertEquals("2S8M", readInfo.ConsensusRead.getCigarString());
        assertEquals("TTCGATAAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(13, readInfo.ConsensusRead.getAlignmentStart());
    }

    @Test
    public void testHardClippedReads()
    {
        String readBases = "CTTCGATAATGGCCG";
        SAMRecord read1 = createSamRecord(nextReadId(), 13, readBases, "3H8M7S", false);
        SAMRecord read2 = createSamRecord(nextReadId(), 15, readBases, "5S10M5H", false);
        SAMRecord read3 = createSamRecord(nextReadId(), 12, readBases, "2H13M2S", false);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2, read3), UMI_ID_1);

        assertEquals("2H13M2S", readInfo.ConsensusRead.getCigarString());
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals(12, readInfo.ConsensusRead.getAlignmentStart());

        readBases = "CTTCGATAATGGCCGGGCCG";
        read1 = createSamRecord(nextReadId(), 10, readBases, "5H10M10D10M10H", false);
        read2 = createSamRecord(nextReadId(), 10, readBases, "5H10M10D10M10H", false);
        read3 = createSamRecord(nextReadId(), 15, readBases, "10H5M10D10M10H", false);

        readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2, read3), UMI_ID_1);

        assertEquals("5H10M10D10M10H", readInfo.ConsensusRead.getCigarString());
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());
    }

    @Test
    public void testMismatchedReadsInsert()
    {
        SAMRecord read1 = createSamRecord(nextReadId(), 10, "CTTCGATAAAAT", "7M1I4M", true);
        SAMRecord read2 = createSamRecord(nextReadId(), 13, "TTCGATAAAT", "2S8M", true);

        // they are marked as same read cause read is reversed and the end are the same
        assertEquals(read1.getAlignmentEnd(), read2.getAlignmentEnd());

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2), UMI_ID_1);

        assertEquals("7M1I4M", readInfo.ConsensusRead.getCigarString());
        assertEquals("CTTCGATAAAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(10, readInfo.ConsensusRead.getAlignmentStart());

        // add one more read so the first C will get chopped
        SAMRecord read3 = createSamRecord(nextReadId(), 13, "TTCGATAAAT", "2S8M", true);

        readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2, read3), UMI_ID_1);

        assertEquals("2S8M", readInfo.ConsensusRead.getCigarString());
        assertEquals("TTCGATAAAT", readInfo.ConsensusRead.getReadString());
        assertEquals(13, readInfo.ConsensusRead.getAlignmentStart());
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

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, List.of(read1, read2, read3, read4), UMI_ID_1);

        assertEquals("30S30M", readInfo.ConsensusRead.getCigarString()); // was 30M2D31M when favouring non-SCs
        assertEquals(40, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(69, readInfo.ConsensusRead.getAlignmentEnd());
        String consensusBases = REF_BASES.substring(40, 100);
        assertEquals(consensusBases, readInfo.ConsensusRead.getReadString());
    }

    private String nextReadId() { return nextUmiReadId(UMI_ID_1, mReadIdGen); }

    private static SAMRecord createSamRecord(
            final String readId, int readStart, final String readBases, final String cigar, boolean isReversed)
    {
        return SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, 5000, isReversed, false, null);
    }

}
