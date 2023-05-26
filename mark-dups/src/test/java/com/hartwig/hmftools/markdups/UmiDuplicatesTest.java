package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.alignmentsToSamTag;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.markdups.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.umi.PositionFragmentCounts;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UmiDuplicatesTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final RecordWriter mWriter;

    private final ChrBaseRegion mChrBaseRegion;
    private final ChromosomeReader mChrReaderUMIs;
    private final ChromosomeReader mChrReaderDuplexUMIs;

    public UmiDuplicatesTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();

        MarkDupsConfig umiConfig = new MarkDupsConfig(1000, 1000, mRefGenome, true, false, false);
        mWriter = new RecordWriter(umiConfig);
        mWriter.setCacheReads();

        mChrBaseRegion = new ChrBaseRegion(CHR_1, 1, 100000);

        mChrReaderUMIs = new ChromosomeReader(mChrBaseRegion, umiConfig, mWriter, new PartitionDataStore(umiConfig));

        MarkDupsConfig duplexUmiConfig = new MarkDupsConfig(1000, 1000, mRefGenome, true, true, false);
        mChrReaderDuplexUMIs = new ChromosomeReader(mChrBaseRegion, duplexUmiConfig, mWriter, new PartitionDataStore(duplexUmiConfig));
    }

    @Test
    public void testUmiGroup()
    {
        // 2 primaries in a UMI group, followed by their mates in the same partition and then their supps in a different partition
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.duplicateGroupMap().size());
        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate2);

        mChrReaderUMIs.processRead(mate1);
        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

        mChrReaderUMIs.processRead(mate2);
        assertEquals(4, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(supp1);
        assertEquals(4, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

        mChrReaderUMIs.processRead(supp2);
        mChrReaderUMIs.onChromosomeComplete();
        assertEquals(6, mWriter.recordWriteCount());
        assertEquals(3, mWriter.recordWriteCountConsensus());
        assertTrue(partitionData.duplicateGroupMap().isEmpty());
    }

    @Test
    public void testUmiGroupMultipleSupplementaries()
    {
        // one primary has multiple supplementaries and the other has a different supp
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null);

        List<SupplementaryReadData> read1Supps = Lists.newArrayList(
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_2, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_3, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setAttribute(SUPPLEMENTARY_ATTRIBUTE, alignmentsToSamTag(read1Supps));
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos + 2, TEST_READ_BASES, "2S98M", CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_2, suppPos + 1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.duplicateGroupMap().size());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos + 2, true,
                false, null);
        setSecondInPair(mate2);

        mChrReaderUMIs.processRead(mate1);
        mChrReaderUMIs.processRead(mate2);

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read1Supp1);
        mChrReaderUMIs.onChromosomeComplete();

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_2, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord read2Supp1 = createSamRecord(
                read2.getReadName(), CHR_2, suppPos + 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos + 2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read1Supp2);
        mChrReaderUMIs.processRead(read2Supp1);
        mChrReaderUMIs.onChromosomeComplete();

        SAMRecord read1Supp3 = createSamRecord(
                read1.getReadName(), CHR_3, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read1Supp3);
        mChrReaderUMIs.onChromosomeComplete();

        assertEquals(7, mWriter.recordWriteCount());
        assertEquals(3, mWriter.recordWriteCountConsensus());
        assertTrue(partitionData.duplicateGroupMap().isEmpty());
        assertEquals(1, partitionData.incompleteFragmentMap().size());

        partitionData.writeRemainingReads(mWriter, mChrReaderUMIs.consensusReads(), false);
        assertEquals(8, mWriter.recordWriteCount());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
    }

    @Test
    public void testUmiGroupInconsistentSupplementaries()
    {
        // supplementaries arriving out of order and mapped to different locations
        String umidId = "TATTAT";

        int supp2Pos1 = 1000;
        int readPos = 2000;
        int supp2Pos2 = 3000;
        int matePos = 4000;
        int supp1Pos1 = 5000;
        int supp1Pos2 = 6000;

        String readId1 = nextReadId(umidId);
        String readId2 = nextReadId(umidId);

        SAMRecord read2Supp1 = createSamRecord(
                readId2, CHR_1, supp2Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read2Supp1);

        SAMRecord read1 = createSamRecord(
                readId1, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp1Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                readId2, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp2Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        SAMRecord read2Supp2 = createSamRecord(
                readId2, CHR_1, supp2Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read2Supp2);

        SAMRecord mate1 = createSamRecord(
                readId1, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp1Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                readId2, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp2Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));
        setSecondInPair(mate2);

        mChrReaderUMIs.processRead(mate1);
        mChrReaderUMIs.processRead(mate2);

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read1Supp1);
        mChrReaderUMIs.processRead(read1Supp2);

        mChrReaderUMIs.onChromosomeComplete();

        assertEquals(8, mWriter.recordWriteCount());
        assertEquals(4, mWriter.recordWriteCountConsensus());

        for(PartitionData partitionData : mChrReaderUMIs.partitionDataStore().partitions())
        {
            assertTrue(partitionData.duplicateGroupMap().isEmpty());
            assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        }
    }

    @Test
    public void testSingleUmiGroup()
    {
        // 2 UMI groups and a single fragment with the forward orientation, 1 fragment and UMI group on the reverse
        // are matched off with each other to form 3 final groups in total
        int readPos = 100;
        int matePos = 200;

        String umi1 = "AAAAA";
        String umi2 = "TTTTT";
        String umi3 = "CCCCC";
        String umi4 = "GGGGG";
        String umi5 = "ACGTC";

        SAMRecord read1 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read1);

        SAMRecord read2 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read2);

        SAMRecord read3 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read3);

        SAMRecord read4 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read4);

        SAMRecord read5 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read5);

        SAMRecord read6 = createSamRecord(
                nextReadId(umi3), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read6);

        // and the reverse orientation reads
        SAMRecord read7 = createSamRecord(
                nextReadId(umi4), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read8 = createSamRecord(
                nextReadId(umi5), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read9 = createSamRecord(
                nextReadId(umi5), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.processRead(read3);
        mChrReaderUMIs.processRead(read4);
        mChrReaderUMIs.processRead(read5);
        mChrReaderUMIs.processRead(read6);
        mChrReaderUMIs.processRead(read7);
        mChrReaderUMIs.processRead(read8);
        mChrReaderUMIs.processRead(read9);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.umiGroups().size());
        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(9, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

        assertEquals(2, mChrReaderUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mChrReaderUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mChrReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
        assertEquals(3, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);
    }

    @Test
    public void testSingleUmiGroup2()
    {
        // 2 separate groups of single opposite orientation fragments, where each can be collapsed
        int readPos = 100;
        int matePos = 200;
        int matePos2 = 300;

        String umi1 = "AAAAA";
        String umi2 = "TTTTT";
        String umi3 = "CCCCC";
        String umi4 = "GGGGG";
        String umi5 = "ACGTA";

        SAMRecord read1 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read2);

        SAMRecord read3 = createSamRecord(
                nextReadId(umi3), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos2, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read3);

        SAMRecord read4 = createSamRecord(
                nextReadId(umi4), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos2, false, false,
                null, true, TEST_READ_CIGAR);

        // and one excess read
        SAMRecord read5 = createSamRecord(
                nextReadId(umi5), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos2, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read5);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.processRead(read3);
        mChrReaderUMIs.processRead(read4);
        mChrReaderUMIs.processRead(read5);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.umiGroups().size());
        assertEquals(5, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

        assertEquals(2, mChrReaderUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mChrReaderUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mChrReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(3, positionFragmentCounts.UniqueCoordCount);
        assertEquals(3, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);
    }

    @Test
    public void testSingleUmiGroup3()
    {
        // 2 single fragments with difference coordinates
        int readPos = 100;
        int matePos = 200;
        int matePos2 = 300;

        String umi1 = "AAAAA";
        String umi2 = "TTTTT";

        SAMRecord read1 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos2, false, false,
                null, true, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(0, mWriter.recordWriteCountConsensus());

        PositionFragmentCounts positionFragmentCounts = mChrReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
        assertEquals(2, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);
    }

    @Test
    public void testSingleUmiGroup4()
    {
        // one UMI group after collapsing
        int readPos = 100;
        int matePos = 200;

        String umi1 = "AAAAA";
        String umi2 = "TTTTT";

        SAMRecord read1 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read3 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read3);

        SAMRecord read4 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read4);

        // and one excess read
        SAMRecord read5 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read5);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.processRead(read3);
        mChrReaderUMIs.processRead(read4);
        mChrReaderUMIs.processRead(read5);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(1, partitionData.umiGroups().size());
        assertEquals(5, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

        assertEquals(1, mChrReaderUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mChrReaderUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mChrReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(1, positionFragmentCounts.UniqueCoordCount);
        assertEquals(1, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);
    }

    @Test
    public void testSingleUmiGroup5()
    {
        // one UMI group without collapsing
        int readPos = 100;
        int matePos = 200;

        String umi1 = "AAAAA";

        SAMRecord read1 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(1, partitionData.umiGroups().size());
        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

        assertEquals(1, mChrReaderUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mChrReaderUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mChrReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(1, positionFragmentCounts.UniqueCoordCount);
        assertEquals(1, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);

        String umi2 = "GGGGG";

        SAMRecord read3 = createSamRecord(
                nextReadId(umi2), CHR_1, readPos + 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos + 1, false, false,
                null, true, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read3);
        mChrReaderUMIs.flushReadPositions();

        assertEquals(2, positionFragmentCounts.Frequency);
    }

    @Test
    public void testDuplexUmiGroup()
    {
        // 2 fragments have opposite fragments but are linked by their duplex UMIs
        String umidId1Part1 = "TATTAT";
        String umidId1Part2 = "GCGGCG";

        String umidId2Part1 = "TTAATT";
        String umidId2Part2 = "GGCCGG";

        String umiId1 = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
        String umiId1Reversed = umidId1Part2 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part1.substring(0, 5) + "G"; // 1 base diff is allowed

        String umiId2 = umidId2Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;
        String umiId2Reversed = umidId2Part2.substring(0, 5) + "A" + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part1;

        int readPos = 100;
        int matePos = 200;

        // 2 single fragments
        SAMRecord read1 = createSamRecord(
                nextReadId(umiId1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umiId1Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read2);

        // a pair and a reversed single
        matePos = 300;

        SAMRecord read3 = createSamRecord(
                nextReadId(umiId2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read4 = createSamRecord(
                nextReadId(umiId2Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read4);

        SAMRecord read5 = createSamRecord(
                nextReadId(umiId2Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read5);

        mChrReaderDuplexUMIs.processRead(read1);
        mChrReaderDuplexUMIs.processRead(read2);
        mChrReaderDuplexUMIs.processRead(read3);
        mChrReaderDuplexUMIs.processRead(read4);
        mChrReaderDuplexUMIs.processRead(read5);
        mChrReaderDuplexUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderDuplexUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(5, partitionData.duplicateGroupMap().size());
        assertEquals(5, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

        assertEquals(2  , mChrReaderDuplexUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mChrReaderDuplexUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mChrReaderDuplexUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
        assertEquals(2, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);
    }

    private String nextReadId(final String umiId)
    {
        return format("%s:%s", mReadIdGen.nextId(), umiId);
    }
}
