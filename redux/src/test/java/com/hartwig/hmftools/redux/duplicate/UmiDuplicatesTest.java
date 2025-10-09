package com.hartwig.hmftools.redux.duplicate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.READ_UNMAPPER_DISABLED;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createPartitionRead;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.redux.ReduxConstants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.redux.consensus.ConsensusReads.formConsensusReadId;
import static com.hartwig.hmftools.redux.duplicate.UmiGroupBuilder.buildUmiGroups;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.PartitionReader;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.TestBamWriter;
import com.hartwig.hmftools.redux.consensus.TemplateReads;

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

        ReduxConfig umiConfig = new ReduxConfig(mRefGenome, true, false, false, READ_UNMAPPER_DISABLED);

        mWriter = new TestBamWriter(umiConfig);

        mPartitionReaderUMIs = createPartitionRead(umiConfig, mWriter);

        ReduxConfig duplexUmiConfig = new ReduxConfig(mRefGenome, true, true, false, READ_UNMAPPER_DISABLED);
        mPartitionReaderDuplexUMIs = createPartitionRead(duplexUmiConfig, mWriter);
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

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        mPartitionReaderUMIs.processRead(read1);
        mPartitionReaderUMIs.processRead(read2);
        mPartitionReaderUMIs.flushReadPositions();

        assertEquals(2, mWriter.nonConsensusWriteCount());
        assertEquals(1, mWriter.consensusWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate2);

        mPartitionReaderUMIs.processRead(mate1);
        assertEquals(2, mWriter.nonConsensusWriteCount());
        assertEquals(1, mWriter.consensusWriteCount());

        mPartitionReaderUMIs.processRead(mate2);

        mPartitionReaderUMIs.flushReadPositions();
        assertEquals(4, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReaderUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        mPartitionReaderUMIs.processRead(supp1);
        assertEquals(4, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

        mPartitionReaderUMIs.processRead(supp2);
        mPartitionReaderUMIs.postProcessRegion();

        assertEquals(6, mWriter.nonConsensusWriteCount());
        assertEquals(3, mWriter.consensusWriteCount());
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

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read1);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read2);

        SAMRecord read3 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read3);

        SAMRecord read4 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read4);

        SAMRecord read5 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read5);

        SAMRecord read6 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi3), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read6);

        // and the reverse orientation reads
        SAMRecord read7 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi4), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read8 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi5), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read9 = SamRecordTestUtils.createSamRecord(
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

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
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
        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId1Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read2);

        // a pair and a reversed single
        matePos = 300;

        SAMRecord read3 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId2), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        SAMRecord read4 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId2Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read4);

        SAMRecord read5 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId2Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);
        setSecondInPair(read5);

        mPartitionReaderDuplexUMIs.processRead(read1);
        mPartitionReaderDuplexUMIs.processRead(read2);
        mPartitionReaderDuplexUMIs.processRead(read3);
        mPartitionReaderDuplexUMIs.processRead(read4);
        mPartitionReaderDuplexUMIs.processRead(read5);
        mPartitionReaderDuplexUMIs.postProcessRegion();

        // PartitionData partitionData = mPartitionReaderDuplexUMIs.partitionDataStore().getOrCreatePartitionData("1_0");
        // assertEquals(5, partitionData.duplicateGroupMap().size());
        assertEquals(5, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

        assertEquals(2, mPartitionReaderDuplexUMIs.statistics().UmiStats.UmiGroups);
        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().UmiStats.PositionFragments.size());

        PositionFragmentCounts positionFragmentCounts = mPartitionReaderDuplexUMIs.statistics().UmiStats.PositionFragments.get(0);
        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
        assertEquals(2, positionFragmentCounts.UniqueFragmentCount);
        assertEquals(1, positionFragmentCounts.Frequency);

        assertEquals(2, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.size());
        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(2).DualStrandCount);
        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(3).DualStrandCount);
    }

    @Test
    public void testAlternativeFirstInPairGroups()
    {
        // 4 primary reads with alternating first in pair designation for the same lower / upper reads
        int readPos = 100;
        int matePos = 1500;
        int suppPos = 2200;
        int suppPos2 = 2400;

        String umidId1Part1 = "TATTAT";
        String umidId1Part2 = "GCGGCG";
        String umiId1 = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
        String umiId1Reversed = umidId1Part2 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part1;

        SupplementaryReadData suppData1 = new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);
        SupplementaryReadData suppData2 = new SupplementaryReadData(CHR_1, suppPos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);

        SAMRecord read1 = createSamRecord(
                nextReadId(umiId1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                false, suppData1, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umiId1Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                false, null, true, TEST_READ_CIGAR);
        setSecondInPair(read2);

        SAMRecord read3 = createSamRecord(
                nextReadId(umiId1), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                false, suppData1, true, TEST_READ_CIGAR);

        SAMRecord read4 = createSamRecord(
                nextReadId(umiId1Reversed), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                false, suppData1, true, TEST_READ_CIGAR);
        setSecondInPair(read4);

        mPartitionReaderDuplexUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        mPartitionReaderDuplexUMIs.processRead(read1);
        mPartitionReaderDuplexUMIs.processRead(read2);
        mPartitionReaderDuplexUMIs.processRead(read3);
        mPartitionReaderDuplexUMIs.processRead(read4);
        mPartitionReaderDuplexUMIs.flushReadPositions();
        mPartitionReaderDuplexUMIs.postProcessRegion();

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, suppData2, false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, suppData2, false, TEST_READ_CIGAR);

        SAMRecord mate3 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate3);

        SAMRecord mate4 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null, false, TEST_READ_CIGAR);

        mPartitionReaderDuplexUMIs.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        mPartitionReaderDuplexUMIs.processRead(mate1);
        mPartitionReaderDuplexUMIs.processRead(mate2);
        mPartitionReaderDuplexUMIs.processRead(mate3);
        mPartitionReaderDuplexUMIs.processRead(mate4);
        mPartitionReaderDuplexUMIs.flushReadPositions();
        mPartitionReaderDuplexUMIs.postProcessRegion();

        // finally the expected supplementaries

        mPartitionReaderDuplexUMIs.setupRegion(new ChrBaseRegion(CHR_1, 2001, 3000));

        SupplementaryReadData suppDataRead1 = new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);
        SupplementaryReadData suppDataRead2 = new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, suppDataRead1, true, TEST_READ_CIGAR);

        // read 2 has no supp

        SAMRecord supp3 = createSamRecord(
                read3.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, suppDataRead1, true, TEST_READ_CIGAR);

        SAMRecord supp4 = createSamRecord(
                read4.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, suppDataRead1, true, TEST_READ_CIGAR);
        setSecondInPair(supp4);

        SAMRecord mateSupp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                true, suppDataRead2, false, TEST_READ_CIGAR);
        setSecondInPair(mateSupp1);

        SAMRecord mateSupp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                true, suppDataRead2, false, TEST_READ_CIGAR);

        mPartitionReaderDuplexUMIs.processRead(supp1);
        mPartitionReaderDuplexUMIs.processRead(supp3);
        mPartitionReaderDuplexUMIs.processRead(supp4);
        mPartitionReaderDuplexUMIs.processRead(mateSupp1);
        mPartitionReaderDuplexUMIs.processRead(mateSupp2);
        mPartitionReaderDuplexUMIs.postProcessRegion();

        assertEquals(13, mWriter.nonConsensusWriteCount());
        assertEquals(4, mWriter.consensusWriteCount());
    }

    @Test
    public void testConsistentReadPairsAfterUmiCollapsing()
    {
        ReduxConfig config = new ReduxConfig(null, true, true, false, READ_UNMAPPER_DISABLED);
        UmiConfig umiConfig = config.UMIs;

        String umiId1 = "GTCACTC" + DEFAULT_DUPLEX_UMI_DELIM + "TAAGATA";
        String umiId2 = "GTCACTC" + DEFAULT_DUPLEX_UMI_DELIM + "TACGATC";

        int matePosition = 200;

        SAMRecord read1 = createSamRecord(
                "READ_001:" + umiId1, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePosition,
                false, false, null, true, TEST_READ_CIGAR);
        SAMRecord mate1 = createSamRecord(
                "READ_001:" + umiId1, CHR_2, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePosition,
                true, false, null, false, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                "READ_002:" + umiId2, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePosition,
                false, false, null, true, TEST_READ_CIGAR);
        SAMRecord mate2 = createSamRecord(
                "READ_002:" + umiId2, CHR_2, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePosition,
                true, false, null, false, TEST_READ_CIGAR);

        DuplicateGroup readGroup = new DuplicateGroup(Lists.newArrayList(read1, read2), FragmentCoords.fromRead(read1, true));
        DuplicateGroup mateGroup = new DuplicateGroup(Lists.newArrayList(mate2, mate1), FragmentCoords.fromRead(mate1, true));

        List<DuplicateGroup> readUmiGroups = buildUmiGroups(readGroup.fragmentCoordinates(), readGroup.reads(), umiConfig);
        List<DuplicateGroup> mateUmiGroups = buildUmiGroups(mateGroup.fragmentCoordinates(), mateGroup.reads(), umiConfig);

        assertEquals(1, readUmiGroups.size());
        assertEquals(1, mateUmiGroups.size());

        DuplicateGroup readUmiGroup = readUmiGroups.get(0);
        DuplicateGroup mateUmiGroup = mateUmiGroups.get(0);

        SAMRecord readTemplate = TemplateReads.selectTemplateRead(readUmiGroup.reads(), readUmiGroup.fragmentCoordinates());
        String readConsensusReadName = formConsensusReadId(readTemplate, readUmiGroup.umi());

        SAMRecord mateTemplate = TemplateReads.selectTemplateRead(mateUmiGroup.reads(), mateUmiGroup.fragmentCoordinates());
        String mateConsensusReadName = formConsensusReadId(mateTemplate, mateUmiGroup.umi());

        assertEquals(readConsensusReadName, mateConsensusReadName);
    }

    @Test
    public void testPolyGCollapsingCreatedOrphanConsensusRead()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(2_000));
        refGenome.ChromosomeLengths.put(CHR_1, 2_000);
        refGenome.RefGenomeMap.put(CHR_2, "A".repeat(2_000));
        refGenome.ChromosomeLengths.put(CHR_2, 2_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        String umiIdPart1 = "A".repeat(7);
        String umiId1Part2 = "G".repeat(7);
        String umiId2Part2 = "C" + "G".repeat(6);

        String readName1 = "READ_001:" + umiIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiId1Part2;
        String readName2 = "READ_002:" + umiIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiId2Part2;

        int readPosition = 1000;
        int matePosition = 1000;

        SAMRecord read1 = createSamRecord(readName1, CHR_1, readPosition, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePosition,
                false, false, null, false, TEST_READ_CIGAR);
        flipFirstInPair(read1);

        SAMRecord mate1 = createSamRecord(readName1, CHR_2, readPosition, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePosition,
                false, false, null, false, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(readName2, CHR_2, readPosition, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePosition,
                false, false, null, false, NO_CIGAR);
        read2.setMateUnmappedFlag(true);

        SAMRecord mate2 = createSamRecord(readName2, CHR_2, readPosition, TEST_READ_BASES, NO_CIGAR, CHR_2, matePosition,
                false, false, null, false, TEST_READ_CIGAR);
        flipFirstInPair(mate2);
        mate2.setReadUnmappedFlag(true);

        List<SAMRecord> reads = Lists.newArrayList(read1, mate1, read2, mate2);
        Map<String, List<SAMRecord>> readsByChromosome = Maps.newHashMap();

        for(SAMRecord read : reads)
        {
            String chromosome = read.getReferenceName();
            readsByChromosome.computeIfAbsent(chromosome, k -> Lists.newArrayList());
            readsByChromosome.get(chromosome).add(read);
        }

        for(List<SAMRecord> chromosomeReads : readsByChromosome.values())
            chromosomeReads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));

        for(String chromosome : readsByChromosome.keySet())
        {
            partitionReader.setupRegion(new ChrBaseRegion(chromosome, 1, 2_000));
            readsByChromosome.get(chromosome).forEach(x -> partitionReader.processRead(x.deepCopy()));
            partitionReader.postProcessRegion();
        }

        assertEquals(4, writer.nonConsensusWriteCount());
        assertEquals(0, writer.consensusWriteCount());

        final List<SAMRecord> writtenRecords = Collections.unmodifiableList(writer.WrittenRecords);
        final BiFunction<String, Boolean, Stream<SAMRecord>>
                readStream = (final String readName, final Boolean isFirstOfPair) -> writtenRecords.stream()
                .filter(x -> (x.getReadName().equals(readName) && x.getFirstOfPairFlag() == isFirstOfPair));

        assertEquals(1L, readStream.apply(readName1, true).count());
        assertEquals(1L, readStream.apply(readName1, false).count());
        assertEquals(1L, readStream.apply(readName2, true).count());
        assertEquals(1L, readStream.apply(readName2, false).count());

        SAMRecord outRead1 = readStream.apply(readName1, false).findAny().orElse(null);
        SAMRecord outMate1 = readStream.apply(readName1, true).findAny().orElse(null);
        SAMRecord outRead2 = readStream.apply(readName2, true).findAny().orElse(null);
        SAMRecord outMate2 = readStream.apply(readName2, false).findAny().orElse(null);

        outRead1.setAttribute(UMI_ATTRIBUTE, null);
        outMate1.setAttribute(UMI_ATTRIBUTE, null);
        outRead2.setAttribute(UMI_ATTRIBUTE, null);
        outMate2.setAttribute(UMI_ATTRIBUTE, null);

        assertEquals(read1, outRead1);
        assertEquals(mate1, outMate1);

        assertTrue(outRead2.getDuplicateReadFlag());
        outRead2.setDuplicateReadFlag(false);
        assertEquals(read2, outRead2);

        assertTrue(outMate2.getDuplicateReadFlag());
        outMate2.setDuplicateReadFlag(false);
        assertEquals(mate2, outMate2);
    }

    @Test
    public void testPolyGCollapsingOnUnmappedReads()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(2_000));
        refGenome.ChromosomeLengths.put(CHR_1, 2_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        String umiIdPart1 = "A".repeat(7);
        String umiId1Part2 = "G".repeat(7);
        String umiId2Part2 = "A".repeat(7);

        String readName1 = "READ_001:" + umiIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiId1Part2;
        String readName2 = "READ_002:" + umiIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiId1Part2;
        String readName3 = "READ_003:" + umiIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiId2Part2;
        String readName4 = "READ_004:" + umiIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiId2Part2;

        SAMRecord read1 = createSamRecord(
                readName1, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false,
                false, null, false, TEST_READ_CIGAR);
        read1.setReadUnmappedFlag(true);

        SAMRecord mate1 = createSamRecord(
                readName1, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, false, NO_CIGAR);
        mate1.setMateUnmappedFlag(true);
        flipFirstInPair(mate1);

        SAMRecord read2 = createSamRecord(
                readName2, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false,
                false, null, false, TEST_READ_CIGAR);
        read2.setReadUnmappedFlag(true);

        SAMRecord mate2 = createSamRecord(
                readName2, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, false, NO_CIGAR);
        mate2.setMateUnmappedFlag(true);
        flipFirstInPair(mate2);

        SAMRecord read3 = createSamRecord(
                readName3, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false,
                false, null, false, TEST_READ_CIGAR);
        read3.setReadUnmappedFlag(true);

        SAMRecord mate3 = createSamRecord(
                readName3, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, false, NO_CIGAR);
        mate3.setMateUnmappedFlag(true);
        flipFirstInPair(mate3);

        SAMRecord read4 = createSamRecord(
                readName4, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false,
                false, null, false, TEST_READ_CIGAR);
        read4.setReadUnmappedFlag(true);

        SAMRecord mate4 = createSamRecord(
                readName1, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, false, NO_CIGAR);
        mate4.setMateUnmappedFlag(true);
        flipFirstInPair(mate4);

        List<SAMRecord> reads = List.of(mate1, mate2, mate3, mate4, read1, read2, read3, read4);
        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 2_000));
        reads.forEach(partitionReader::processRead);
        partitionReader.postProcessRegion();

        assertEquals(8, writer.nonConsensusWriteCount());
        assertEquals(2, writer.consensusWriteCount());
    }

    @Test
    public void testPolyGConsistentCollapsingOfReadsAndTheirUnmappedMates()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(2_000));
        refGenome.ChromosomeLengths.put(CHR_1, 2_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        final PartitionReader partitionReader = createPartitionRead(config, writer);

        String umiId = "A".repeat(7) + DEFAULT_DUPLEX_UMI_DELIM + "G".repeat(7);

        String readName1 = nextReadId(umiId);
        String readName2 = nextReadId(umiId);
        String readName3 = nextReadId(umiId);

        SAMRecord read1 = createSamRecord(
                readName1, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false,
                false, null, false, NO_CIGAR);
        read1.setMateUnmappedFlag(true);

        SAMRecord mate1 = createSamRecord(
                readName1, CHR_1, 500, TEST_READ_BASES, NO_CIGAR, CHR_1, 500, false,
                false, null, false, TEST_READ_CIGAR);
        mate1.setReadUnmappedFlag(true);
        flipFirstInPair(mate1);

        SAMRecord read2 = createSamRecord(
                readName2, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false,
                false, null, false, NO_CIGAR);
        read2.setMateUnmappedFlag(true);

        SAMRecord mate2 = createSamRecord(
                readName2, CHR_1, 500, TEST_READ_BASES, NO_CIGAR, CHR_1, 500, false,
                false, null, false, TEST_READ_CIGAR);
        mate2.setReadUnmappedFlag(true);
        flipFirstInPair(mate2);

        SAMRecord read3 = createSamRecord(
                readName3, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_500, false,
                false, null, false, TEST_READ_CIGAR);
        SAMRecord mate3 = createSamRecord(
                readName3, CHR_1, 1_500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false,
                false, null, false, TEST_READ_CIGAR);
        flipFirstInPair(mate3);

        List<SAMRecord> reads = Lists.newArrayList(read1, mate1, read2, mate2, read3, mate3);
        reads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 2_000));
        reads.forEach(x -> partitionReader.processRead(x.deepCopy()));
        partitionReader.postProcessRegion();

        assertEquals(6, writer.nonConsensusWriteCount());
        assertEquals(0, writer.consensusWriteCount());

        final Function<SAMRecord, SAMRecord> makeReadDuplicate = (final SAMRecord read) ->
        {
            SAMRecord readCopy = read.deepCopy();
            readCopy.setDuplicateReadFlag(true);
            return readCopy;
        };

        final Function<SAMRecord, SAMRecord> removeUmiAttribute = (final SAMRecord read) ->
        {
            SAMRecord readCopy = read.deepCopy();
            readCopy.setAttribute(UMI_ATTRIBUTE, null);
            return readCopy;
        };

        SortedMultiset<String> expectedDuplicateReads = Lists.newArrayList(read1, mate1, read2, mate2).stream()
                .map(makeReadDuplicate)
                .map(SAMRecord::getSAMString)
                .collect(Collectors.toCollection(TreeMultiset::create));
        SortedMultiset<String> actualDuplicateReads = writer.WrittenRecords.stream()
                .filter(SAMRecord::getDuplicateReadFlag)
                .map(removeUmiAttribute)
                .map(SAMRecord::getSAMString)
                .collect(Collectors.toCollection(TreeMultiset::create));

        assertEquals(expectedDuplicateReads, actualDuplicateReads);

        // now create a dup of read3
        writer = new TestBamWriter(config);
        final PartitionReader partitionReader2 = createPartitionRead(config, writer);

        String readName4 = nextReadId(umiId);
        SAMRecord read4 = createSamRecord(
                readName4, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_500, false,
                false, null, false, TEST_READ_CIGAR);
        SAMRecord mate4 = createSamRecord(
                readName4, CHR_1, 1_500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false,
                false, null, false, TEST_READ_CIGAR);
        flipFirstInPair(mate4);

        reads.add(read4);
        reads.add(mate4);
        reads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));

        partitionReader2.setupRegion(new ChrBaseRegion(CHR_1, 1, 2_000));
        reads.forEach(x -> partitionReader2.processRead(x.deepCopy()));
        partitionReader2.postProcessRegion();

        assertEquals(8, writer.nonConsensusWriteCount());
        assertEquals(2, writer.consensusWriteCount());

        expectedDuplicateReads = Lists.newArrayList(read1, mate1, read2, mate2, read3, mate3, read4, mate4).stream()
                .map(makeReadDuplicate)
                .map(SAMRecord::getSAMString)
                .collect(Collectors.toCollection(TreeMultiset::create));
        actualDuplicateReads = writer.WrittenRecords.stream()
                .filter(SAMRecord::getDuplicateReadFlag)
                .map(removeUmiAttribute)
                .map(SAMRecord::getSAMString)
                .collect(Collectors.toCollection(TreeMultiset::create));

        assertEquals(expectedDuplicateReads, actualDuplicateReads);
    }

    protected static String nextReadId(final String umiId)
    {
        return format("%s:%s", READ_ID_GEN.nextId(), umiId);
    }
}
