package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.alignmentsToSamTag;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.redux.common.PartitionData;
import com.hartwig.hmftools.redux.umi.PositionFragmentCounts;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UmiDuplicatesTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final TestBamWriter mWriter;

    private final PartitionReader mPartitionReaderUMIs;
    private final PartitionReader mPartitionReaderDuplexUMIs;

    public UmiDuplicatesTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES_REPEAT_40);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES_REPEAT_40.length());

        ReduxConfig umiConfig = new ReduxConfig(1000, 1000, mRefGenome, true, false, false);

        mWriter = new TestBamWriter(umiConfig);

        mPartitionReaderUMIs = new PartitionReader(umiConfig, null, mWriter, new PartitionDataStore(umiConfig));

        ReduxConfig duplexUmiConfig = new ReduxConfig(1000, 1000, mRefGenome, true, true, false);
        mPartitionReaderDuplexUMIs = new PartitionReader(duplexUmiConfig, null, mWriter, new PartitionDataStore(duplexUmiConfig));
    }

    @Test
    public void testUmiGroup()
    {
        // 2 primaries in a UMI group, followed by their mates in the same partition and then their supps in a different partition
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int suppPos = 1800;

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        SAMRecord read1 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        mPartitionReaderUMIs.processRead(read1);
        mPartitionReaderUMIs.processRead(read2);
        mPartitionReaderUMIs.flushReadPositions();

        PartitionData partitionData = mPartitionReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.duplicateGroupMap().size());
        assertEquals(2, mWriter.nonConsensusWriteCount());
        assertEquals(1, mWriter.consensusWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate2);

        mPartitionReaderUMIs.processRead(mate1);
        assertEquals(2, mWriter.nonConsensusWriteCount());
        assertEquals(1, mWriter.consensusWriteCount());

        mPartitionReaderUMIs.processRead(mate2);
        assertEquals(4, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        mPartitionReaderUMIs.processRead(supp1);
        assertEquals(4, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

        mPartitionReaderUMIs.processRead(supp2);
        mPartitionReaderUMIs.postProcessRegion();

        assertEquals(6, mWriter.nonConsensusWriteCount());
        assertEquals(3, mWriter.consensusWriteCount());
        assertTrue(partitionData.duplicateGroupMap().isEmpty());
    }

    @Test
    public void testUmiGroupMultipleSupplementaries()
    {
        // one primary has multiple supplementaries and the other has a different supp
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int suppPos = 1800;

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

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        mPartitionReaderUMIs.processRead(read1);
        mPartitionReaderUMIs.processRead(read2);
        mPartitionReaderUMIs.flushReadPositions();

        PartitionData partitionData = mPartitionReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.duplicateGroupMap().size());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos + 2, true,
                false, null);
        setSecondInPair(mate2);

        mPartitionReaderUMIs.processRead(mate1);
        mPartitionReaderUMIs.processRead(mate2);

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReaderUMIs.processRead(read1Supp1);
        mPartitionReaderUMIs.postProcessRegion();

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_2, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord read2Supp1 = createSamRecord(
                read2.getReadName(), CHR_2, suppPos + 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos + 2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_2, 1001, 2000));

        mPartitionReaderUMIs.processRead(read1Supp2);
        mPartitionReaderUMIs.processRead(read2Supp1);
        mPartitionReaderUMIs.postProcessRegion();

        SAMRecord read1Supp3 = createSamRecord(
                read1.getReadName(), CHR_3, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_3, 1001, 2000));

        mPartitionReaderUMIs.processRead(read1Supp3);
        mPartitionReaderUMIs.postProcessRegion();

        assertEquals(7, mWriter.nonConsensusWriteCount());
        assertEquals(3, mWriter.consensusWriteCount());
        assertTrue(partitionData.duplicateGroupMap().isEmpty());
        assertEquals(1, partitionData.incompleteFragmentMap().size());

        partitionData.writeRemainingReads(mWriter, mPartitionReaderUMIs.consensusReads(), false);
        assertEquals(8, mWriter.nonConsensusWriteCount());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
    }

    @Test
    public void testUmiGroupInconsistentSupplementaries()
    {
        // supplementaries arriving out of order and mapped to different locations
        String umidId = "TATTAT";

        int supp2Pos1 = 950;
        int readPos = 1950;
        int supp2Pos2 = 2950;
        int matePos = 3500;
        int supp1Pos1 = 4500;
        int supp1Pos2 = 5500;

        String readId1 = nextReadId(umidId);
        String readId2 = nextReadId(umidId);

        SAMRecord read2Supp1 = createSamRecord(
                readId2, CHR_1, supp2Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1, 999));
        mPartitionReaderUMIs.processRead(read2Supp1);
        mPartitionReaderUMIs.postProcessRegion();

        SAMRecord read1 = createSamRecord(
                readId1, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp1Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                readId2, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp2Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        int regionStart = 1000;
        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReaderUMIs.processRead(read1);
        mPartitionReaderUMIs.processRead(read2);

        // mChrReaderUMIs.flushReadPositions();
        mPartitionReaderUMIs.postProcessRegion();

        SAMRecord read2Supp2 = createSamRecord(
                readId2, CHR_1, supp2Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        regionStart += 1000;
        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReaderUMIs.processRead(read2Supp2);
        mPartitionReaderUMIs.postProcessRegion();

        SAMRecord mate1 = createSamRecord(
                readId1, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp1Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                readId2, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp2Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));
        setSecondInPair(mate2);

        regionStart += 1000;
        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReaderUMIs.processRead(mate1);
        mPartitionReaderUMIs.processRead(mate2);
        mPartitionReaderUMIs.postProcessRegion();

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        // note last 2 regions are out of order
        regionStart = 5000;
        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReaderUMIs.processRead(read1Supp1);
        mPartitionReaderUMIs.postProcessRegion();

        regionStart = 4000;
        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReaderUMIs.processRead(read1Supp2);
        mPartitionReaderUMIs.postProcessRegion();

        assertEquals(8, mWriter.nonConsensusWriteCount());
        assertEquals(4, mWriter.consensusWriteCount());

        for(PartitionData partitionData : mPartitionReaderUMIs.partitionDataStore().partitions())
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

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

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

        mPartitionReaderUMIs.processRead(read1);
        mPartitionReaderUMIs.processRead(read2);
        mPartitionReaderUMIs.processRead(read3);
        mPartitionReaderUMIs.processRead(read4);
        mPartitionReaderUMIs.processRead(read5);
        mPartitionReaderUMIs.processRead(read6);
        mPartitionReaderUMIs.processRead(read7);
        mPartitionReaderUMIs.processRead(read8);
        mPartitionReaderUMIs.processRead(read9);
        mPartitionReaderUMIs.postProcessRegion();

        PartitionData partitionData = mPartitionReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(3, partitionData.umiGroups().size());
        assertEquals(2, partitionData.resolvedFragmentStateMap().size());
        assertEquals(9, mWriter.nonConsensusWriteCount());
        assertEquals(3, mWriter.consensusWriteCount());

        assertEquals(3, mPartitionReaderUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mPartitionReaderUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mPartitionReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
        assertEquals(5, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);
    }

    @Test
    public void testSingleUmiGroup2()
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

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        mPartitionReaderUMIs.processRead(read1);
        mPartitionReaderUMIs.processRead(read2);
        mPartitionReaderUMIs.postProcessRegion();

        assertEquals(2, mWriter.nonConsensusWriteCount());
        assertEquals(0, mWriter.consensusWriteCount());

        PositionFragmentCounts positionFragmentCounts = mPartitionReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
        assertEquals(2, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);
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

        mPartitionReaderDuplexUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

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

        mPartitionReaderDuplexUMIs.processRead(read1);
        mPartitionReaderDuplexUMIs.processRead(read2);
        mPartitionReaderDuplexUMIs.processRead(read3);
        mPartitionReaderDuplexUMIs.processRead(read4);
        mPartitionReaderDuplexUMIs.processRead(read5);
        mPartitionReaderDuplexUMIs.postProcessRegion();

        PartitionData partitionData = mPartitionReaderDuplexUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(5, partitionData.duplicateGroupMap().size());
        assertEquals(5, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

        assertEquals(2, mPartitionReaderDuplexUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mPartitionReaderDuplexUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
        assertEquals(2, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);

        assertEquals(2, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.size());
        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(2).DualStrandFrequency);
        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(3).DualStrandFrequency);
    }

    private String nextReadId(final String umiId)
    {
        return format("%s:%s", mReadIdGen.nextId(), umiId);
    }
}
