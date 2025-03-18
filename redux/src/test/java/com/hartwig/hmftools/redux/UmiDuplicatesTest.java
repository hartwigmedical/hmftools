package com.hartwig.hmftools.redux;

import static java.lang.String.format;
import static java.util.regex.Matcher.quoteReplacement;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_14;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.SAM_DICTIONARY_V37;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_UNMAPPER_DISABLED;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createPartitionRead;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

import java.util.List;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.jupiter.api.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

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

        SAMRecord read1 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
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

        assertEquals(9, mWriter.nonConsensusWriteCount());
        assertEquals(3, mWriter.consensusWriteCount());

//        assertEquals(3, mPartitionReaderUMIs.statistics().UmiStats.UmiGroups);
//        assertEquals(1, mPartitionReaderUMIs.statistics().UmiStats.PositionFragments.size());
//
//        PositionFragmentCounts positionFragmentCounts = mPartitionReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
//        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
//        assertEquals(5, positionFragmentCounts.UniqueFragmentCount);
//        assertEquals(1, positionFragmentCounts.Frequency);
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

//        PositionFragmentCounts positionFragmentCounts = mPartitionReaderUMIs.statistics().UmiStats.PositionFragments.get(0);
//        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
//        assertEquals(2, positionFragmentCounts.UniqueFragmentCount);
//        assertEquals(1, positionFragmentCounts.Frequency);
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

        // PartitionData partitionData = mPartitionReaderDuplexUMIs.partitionDataStore().getOrCreatePartitionData("1_0");
        // assertEquals(5, partitionData.duplicateGroupMap().size());
        assertEquals(5, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

//        assertEquals(2, mPartitionReaderDuplexUMIs.statistics().UmiStats.UmiGroups);
//        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().UmiStats.PositionFragments.size());

//        PositionFragmentCounts positionFragmentCounts = mPartitionReaderDuplexUMIs.statistics().UmiStats.PositionFragments.get(0);
//        assertEquals(2, positionFragmentCounts.UniqueCoordCount);
//        assertEquals(2, positionFragmentCounts.UniqueFragmentCount);
//        assertEquals(1, positionFragmentCounts.Frequency);

//        assertEquals(2, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.size());
//        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(2).DualStrandFrequency);
//        assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(3).DualStrandFrequency);
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
    public void testIlluminaJitterUmiGroupCollapse()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, false, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

        String umi = "AAAAA";

        SAMRecord read1 = createSamRecord(
                nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umi), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
        assertEquals(1, writer.consensusWriteCount());
    }

    @Test
    public void testIlluminaMateJitterUmiGroupCollapse()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, false, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

        String umi = "AAAAA";

        SAMRecord read1 = createSamRecord(
                nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false,
                false, null, true, TEST_READ_CIGAR);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
        assertEquals(1, writer.consensusWriteCount());
    }

    @Test
    public void testIlluminaJitterDuplexUmiGroupCollapse()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

        String umidIdPart1 = "TATTAT";
        String umidIdPart2 = "GCGGCG";
        String umiId = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart2;
        String umiIdReversed = umidIdPart2 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart1;

        SAMRecord read1 = createSamRecord(
                nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 = createSamRecord(nextReadId(umiIdReversed), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE,
                TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
        read2.setFirstOfPairFlag(false);
        read2.setSecondOfPairFlag(true);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
        assertEquals(1, writer.consensusWriteCount());
    }

    @Test
    public void testIlluminaMateJitterDuplexUmiGroupCollapse()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

        String umidIdPart1 = "TATTAT";
        String umidIdPart2 = "GCGGCG";
        String umiId = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart2;
        String umiIdReversed = umidIdPart2 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart1;

        SAMRecord read1 = createSamRecord(
                nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 = createSamRecord(nextReadId(umiIdReversed), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false, false, null, true, TEST_READ_CIGAR);
        read2.setFirstOfPairFlag(false);
        read2.setSecondOfPairFlag(true);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
        assertEquals(1, writer.consensusWriteCount());
    }

    // TODO: rename
    @Test
    public void testFoo1()
    {
        String chromosome = CHR_1;
        int minReadStart = 100;
        int readStart1 = minReadStart;
        int readStart2 = minReadStart + 8;
        int mateStart = minReadStart + 381;

        int readLength = 143;
        String readBases = "A".repeat(readLength);
        String cigar = readLength + "M";
        int chromosomeLength = mateStart + readLength + 100;

        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(chromosome, "A".repeat(chromosomeLength));
        refGenome.ChromosomeLengths.put(chromosome, chromosomeLength);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(chromosome, 1, chromosomeLength));

        BiFunction<String, String, String> nextDuplexReadId = (umiPart1, umiPart2) -> umiPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiPart2;

        SAMRecord read1 = createSamRecord(nextDuplexReadId.apply("AGCATCT", "CAGTACT"), chromosome, readStart1, readBases, cigar, chromosome, mateStart, false, false, null, true, cigar);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 = createSamRecord(nextDuplexReadId.apply("GTGTCAT", "TGTCGTA"), chromosome, readStart2, readBases, cigar, chromosome, mateStart, false, false, null, true, cigar);
        read2.setFirstOfPairFlag(false);
        read2.setSecondOfPairFlag(true);

        SAMRecord read3 = createSamRecord(nextDuplexReadId.apply("GTGTCAT", "TGTCGTA"), chromosome, readStart2, readBases, cigar, chromosome, mateStart, false, false, null, true, cigar);
        read3.setFirstOfPairFlag(false);
        read3.setSecondOfPairFlag(true);

        SAMRecord read4 = createSamRecord(nextDuplexReadId.apply("TGTCGTA", "GTGTCAT"), chromosome, readStart2, readBases, cigar, chromosome, mateStart, false, false, null, true, cigar);
        read4.setFirstOfPairFlag(true);
        read4.setSecondOfPairFlag(false);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.processRead(read3);
        partitionReader.processRead(read4);
        partitionReader.postProcessRegion();

        assertEquals(4, writer.nonConsensusWriteCount());
        assertEquals(1, writer.consensusWriteCount());
    }

    // TODO: remove
//    @Test
//    public void testFoo2()
//    {
//        String chromosome = CHR_1;
//        int readStart = 100;
//        int mateStart1 = 74;
//        int mateStart2 = 65;
//
//        int readLength = 143;
//        String readBases = "A".repeat(readLength);
//        String cigar = readLength + "M";
//        String mateCigar1 = "30S40M3I70M";
//        String mateCigar2 = "49M21I73M";
//        int chromosomeLength = readStart + readLength + 100;
//
//        MockRefGenome refGenome = new MockRefGenome(true);
//        refGenome.RefGenomeMap.put(chromosome, "A".repeat(chromosomeLength));
//        refGenome.ChromosomeLengths.put(chromosome, chromosomeLength);
//
//        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
//        TestBamWriter writer = new TestBamWriter(config);
//        PartitionReader partitionReader = createPartitionRead(config, writer);
//
//        partitionReader.setupRegion(new ChrBaseRegion(chromosome, 1, chromosomeLength));
//
//        BiFunction<String, String, String> nextDuplexReadId = (umiPart1, umiPart2) -> umiPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiPart2;
//
//        SAMRecord read1 = createSamRecord(nextDuplexReadId.apply("CACTGTC", "TTGGCTC"), chromosome, readStart, readBases, cigar, chromosome, mateStart1, true, false, null, false, mateCigar1);
//        read1.setFirstOfPairFlag(false);
//        read1.setSecondOfPairFlag(true);
//
//        SAMRecord read2 = createSamRecord(nextDuplexReadId.apply("TTGGCTC", "CACTGTC"), chromosome, readStart, readBases, cigar, chromosome, mateStart2, true, false, null, false, mateCigar2);
//        read2.setFirstOfPairFlag(true);
//        read2.setSecondOfPairFlag(false);
//
//        partitionReader.processRead(read1);
//        partitionReader.processRead(read2);
//        partitionReader.postProcessRegion();
//
//        // TODO: keep this around
////        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
////        SAMFileHeader samFileHeader = recordBuilder.getHeader();
////        samFileHeader.setSequenceDictionary(SAM_DICTIONARY_V37);
////        SAMLineParser samLineParser = new SAMLineParser(samFileHeader);
//
//        assertEquals(2, writer.nonConsensusWriteCount());
//        // TODO:
////        assertEquals(1, writer.consensusWriteCount());
//
//        // TODO: remove
//        fail();
//    }

    // TODO: rename
    @Test
//    public void testFoo2()
//    {
//        final List<String> samStrings = Lists.newArrayList(
//                "A01524:289:HFJLTDRX3:1:2101:27733:13322:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCATGAGTCTTATGTCTGATACCAT\tFFFFFF:FF:FFFFFFFFFFFFF:FF:F:FFFF,F:FF:FFFFFFFF,FFFFFFFFFFFFFF:,FFFFFFFFFF:FFFFFF,FFFFFFFFF:F,F:FFFFFF:::FF:FFF:,FFFFFFF,F:FF:,FFFFFFFF,FFFF,FF\tXA:Z:14,+20077871,143M,2;14,-19495446,143M,2;2,-131955600,143M,4;17,+29541306,55M1D88M,6;\tMC:Z:143M\tMD:Z:120C22\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:1:2104:27317:34710:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495446,143M,1;14,+20077871,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2124:10556:20838:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2161:21377:10426:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFF:FFFFFFFFFFFFFF:FFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF:F,FFFFFFFFFFFFF,FFFFFFFF:FF,FFF:FFFFFFFFFFFFFF:,FF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2216:10393:36370:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGACTTATGTCTGATACCAT\tFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFF:F:FFFFFF,:F:FFF:FFFF::F:F,FF,FFFF:FFFFF:FFFFFFF:,:FFFFFFF,,FF,F:F:FFF:F,:F:FFF::F\tXA:Z:14,-19495446,143M,2;14,+20077871,143M,2;\tMC:Z:143M\tMD:Z:125T17\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:1:2229:31656:24721:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFFFFFF:FFF:FFFFFF:FFFFFFFF:FFFFFFFF,FFFFFFFFFFFFFFFFF:FFFFFFFFFFF:F:FFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:,FFFFFFFFFF,FFFFFFF,FFF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2236:29957:8218:TCGTGTA_TGTGCCT\t99\t22\t16345752\t45\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGTTTGGATAGCTATTATACTGAGTCTTATTTCTGATACCAT\tF:FFFFF:FFFF:FFF:F:FFFF,::::,FFFF,FFFF,F:FFFFF:FFFFFFFFFFFFFFFFF,:FF:FFFFFF,F:,FF:FF:F,F:F,FF,FF:FFFFFF,FF,,,FF,FFFFFF,,F,FFFFFF:FF,FFFFFFF,FF,\tXA:Z:14,+20077871,143M,4;14,-19495446,143M,4;\tMC:Z:143M\tMD:Z:103G15C11G11\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:3\tAS:i:128\tXS:i:123",
//                "A01524:289:HFJLTDRX3:1:2238:5023:7623:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFF,FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FF:FFFFFFFF,FFFFFFFFFFFFFF:FFF:FF:FFFF:FFFFF,FFFFFFFFFF:F,,FFFFFFF:FFFF:FFF,FFFFFFF,FFFFFFFFFF:FFFFFFFFF,FF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2244:18810:19132:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495446,143M,1;14,+20077871,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2253:3558:12633:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF,FFF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2255:11966:21245:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:F:FFFFFFFFFFFFFFF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2255:12590:21355:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FF:FFFFFFFFFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2113:14018:28119:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tF::F:FFF,FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFF:FF:FFFFFFFFFF:FFF:FF,:FFFFFFFFFFF,FFFFF,FFFFFFFFF,F:FF,FFFFFF:FFFFF:FFFFFFFF::FFF:,FF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2146:4300:35149:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495446,143M,1;14,+20077871,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2223:7708:8891:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFF:FF,FFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFF:FFFF,FFFFFFFF:FFFFFFFFFF:FFF:FF:F,FF\tXA:Z:14,-19495446,143M,1;14,+20077871,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2237:7166:13714:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFF:FF,FF:FF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFF:FFF:FFFFFFFF:FFFFFF:FF,FFFFFFFF,:FFF,FFFF:FFFFFFF,F:FFFFFFFFFFFFF:FF,FFF\tXA:Z:14,-19495446,143M,1;14,+20077871,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2264:12274:20400:TCGTGTA_TGTGCCT\t99\t22\t16345752\t46\t143M\t=\t16345984\t375\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495446,143M,1;14,+20077871,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCGTGTA_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2274:4860:33739:TGTGCCT_TGTGCCT\t163\t22\t16345752\t46\t143M\t=\t16345918\t309\tCTGACCTTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCAT\tFFFFFFFF:FFFFFFFF,FFF:FFFFF:FFFFFFFFFFFFFFFF,F,FFFFFFFFFFFFFF,FFFFFF,FFFFFFFFFFFF:,F:FF,FFFFFFFFFF:FFF:FFFFFFFFF:FFFFFF,FFFFFFFFFFF,FFFFFFFFFFF\tXA:Z:14,+20077871,143M,1;14,-19495446,143M,1;2,-131955600,143M,3;17,+29541306,55M1D88M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:TGTGCCT_TGTGCCT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2221:19036:15890:CGATATA_CACTGTT\t161\t22\t16345759\t6\t143M\t14\t19495018\t0\tTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTG\tFFFFFFFFFFFFF,F:,FF::F:FFFFFFFF:FFFFF:FF:FFFFFF,FF:FFFF:FF:FFFFFFFFFFFFFFFFF:FFF:F:FFFF:F,FFFFF:::FFFFFFF:FFFFFFFFF:F,F:FFFF:FFFFFFF:FF:,FFFFF,\tXA:Z:14,+20077878,143M,1;14,-19495439,143M,1;2,-131955593,143M,3;17,+29541313,48M1D95M,5;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2161:31286:17691:CGATATA_CACTGTT\t129\t22\t16345759\t6\t143M\t=\t16346180\t422\tTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTG\tFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFF:FFFFFFFFFFFFFFFFFFFFFFF,F,:FFFFF:FFFFFFFFFFFFFF:FFFFFFFFFFFFFF:F:FFFFFFF\tXA:Z:14,-19495439,143M,1;14,+20077878,143M,1;2,-131955593,143M,3;17,+29541313,48M1D95M,5;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2245:13078:32941:CGATATA_CACTGTT\t129\t22\t16345759\t6\t143M\t14\t20078299\t0\tTATGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTG\tFFFFFFFF:F:F:FF:FFFFF:FFFFFF:FFFFFFF:FFFFFFFFFF:FFFFFFF:F::FFFF:FF:FFFF:FFFFFFFFF::FFFFFF:F,FF,FFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF,FFF:FFFFF:\tXA:Z:14,+20077878,143M,1;14,-19495439,143M,1;2,-131955593,143M,3;17,+29541313,48M1D95M,5;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2116:7735:16141:TTGGCTA_AACACAT\t163\t22\t16345762\t6\t143M\t=\t16345906\t287\tGCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATGAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTT\t,FF:FF:FFF:F,::FF:,FF:F,F:F:FF:FFF,FFFF,FFFFFFFFFFFFF,F:FFFFFFFFF:FFF,:,FFFFF:FFFF::FF,FFFFFF,FF,FFFFF,FFFFFFFF:::F:,F:FF,FFFF:FF:,FFFFFFFFFFFF\tXA:Z:14,+20077881,143M,2;14,-19495436,143M,2;2,-131955590,143M,4;17,+29541316,45M1D98M,6;\tMC:Z:143M\tMD:Z:77A65\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:1:2123:21323:30843:TGTGCCT_CAGACGT\t83\t22\t16345763\t46\t143M\t=\t16345665\t-241\tCTTACTATTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTT\tFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-20077882,143M,1;14,+19495435,143M,1;17,-29541317,44M1D99M,5;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2178:30535:8688:GCCATTT_TCCTATG\t99\t22\t16345771\t40\t143M\t=\t16345931\t303\tTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAG\tFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF,FFFFF,FFFFFFFFFF\tXA:Z:14,+20077890,143M,0;14,-19495427,143M,0;2,-131955581,143M,4;17,+29541325,36M1D107M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GCCATTT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2109:9579:12383:GCCATTT_TCCTATG\t99\t22\t16345771\t40\t143M\t=\t16345931\t303\tTGAGTGTTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAG\tFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF\tXA:Z:14,-19495427,143M,0;14,+20077890,143M,0;2,-131955581,143M,4;17,+29541325,36M1D107M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GCCATTT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2162:3848:35744:AGTACAT_AACACAT\t163\t22\t16345778\t25\t143M\t=\t16345903\t268\tTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTA\tFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF:F\tXA:Z:14,+20077897,143M,0;14,-19495420,143M,0;2,-131955574,143M,5;17,+29541332,29M1D114M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:AACACAT_AGTACAT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2149:31430:28902:AACACAT_AGTACAT\t99\t22\t16345778\t25\t143M\t=\t16345903\t268\tTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTA\t:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077897,143M,0;14,-19495420,143M,0;2,-131955574,143M,5;17,+29541332,29M1D114M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:AACACAT_AGTACAT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2175:15664:1752:AACACAT_AGTACAT\t99\t22\t16345778\t25\t143M\t=\t16345903\t268\tTTCTACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFF,FFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,:FFFFFF:FFFFFFFF:FFFFFF\tXA:Z:14,+20077897,143M,0;14,-19495420,143M,0;2,-131955574,143M,5;17,+29541332,29M1D114M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:AACACAT_AGTACAT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2213:17372:20995:TACGATA_TCCTATT\t163\t22\t16345782\t40\t143M\t=\t16345950\t311\tACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATT\tFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF:FFF:F:FFF:FFFFFFFFFFFFFFF,FFFF,FFFFF:,FFF:FF::,FFFF:FFFFFF\tXA:Z:14,+20077901,143M,0;14,-19495416,143M,0;2,-131955570,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCCTATT_TACGATA\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2236:11261:14074:TACGATA_TCCTATT\t163\t22\t16345782\t40\t143M\t=\t16345950\t311\tACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATT\tFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF\tXA:Z:14,+20077901,143M,0;14,-19495416,143M,0;2,-131955570,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCCTATT_TACGATA\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2109:22616:8719:TCCTATT_TACGATA\t99\t22\t16345782\t40\t143M\t=\t16345950\t311\tACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATT\tFFFFFF:F:F:FFF::F::FFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFF:FFF::FFFFFFFFFFF::F:FF:FFFFF,:FFF,F:FF:FFFFFFFFF:FFFFFFF:,:FFF:F:FFFFFFFFF:FFFFF:F:FFF\tXA:Z:14,+20077901,143M,0;14,-19495416,143M,0;2,-131955570,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCCTATT_TACGATA\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2171:21658:1391:TACGATA_TCCTATT\t163\t22\t16345782\t40\t143M\t=\t16345950\t311\tACTTATACCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATT\tFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFF:FFFFFF:FF:FF:FFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF:FFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFF:\tXA:Z:14,-19495416,143M,0;14,+20077901,143M,0;2,-131955570,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:TCCTATT_TACGATA\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2166:6433:36996:CGGTTGT_TCCTATG\t99\t22\t16345789\t40\t143M\t=\t16345963\t317\tCCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAA\tFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF,FFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077908,143M,0;14,-19495409,143M,0;2,-131955563,143M,5;17,+29541343,18M1D125M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2213:31702:25645:CGGTTGT_TCCTATG\t99\t22\t16345789\t40\t143M\t=\t16345963\t317\tCCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAA\tFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFF,FF\tXA:Z:14,+20077908,143M,0;14,-19495409,143M,0;2,-131955563,143M,5;17,+29541343,18M1D125M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2230:22354:29528:CGGTTGT_TCCTATG\t99\t22\t16345789\t40\t143M\t=\t16345963\t317\tCCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAA\tFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFF:FFFFFFFF:FFFFFFFFFFFFF:FFFF,:F\tXA:Z:14,+20077908,143M,0;14,-19495409,143M,0;2,-131955563,143M,5;17,+29541343,18M1D125M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2235:16161:27226:TCCTATG_CGGTTGT\t163\t22\t16345789\t40\t143M\t=\t16345963\t317\tCCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF,:F:,:F:FFF:FF,::F\tXA:Z:14,-19495409,143M,0;14,+20077908,143M,0;2,-131955563,143M,5;17,+29541343,18M1D125M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2125:17219:24361:CGGTTGT_TCCTATG\t99\t22\t16345789\t40\t143M\t=\t16345963\t317\tCCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFF\tXA:Z:14,-19495409,143M,0;14,+20077908,143M,0;2,-131955563,143M,5;17,+29541343,18M1D125M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2153:12979:15671:CGGTTGT_TCCTATG\t99\t22\t16345789\t40\t143M\t=\t16345963\t317\tCCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495409,143M,0;14,+20077908,143M,0;2,-131955563,143M,5;17,+29541343,18M1D125M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2153:13006:15718:CGGTTGT_TCCTATG\t99\t22\t16345789\t40\t143M\t=\t16345963\t317\tCCACACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAA\tFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF:FFFFFFFFFFFFFFFFF\tXA:Z:14,+20077908,143M,0;14,-19495409,143M,0;2,-131955563,143M,5;17,+29541343,18M1D125M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCCTATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2151:24234:22200:GTGTCAT_GCTAACT\t163\t22\t16345792\t40\t143M\t=\t16346018\t369\tCACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FF:FFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:F\tXA:Z:14,+20077911,143M,0;14,-19495406,143M,0;2,-131955560,143M,5;17,+29541346,15M1D128M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGTCAT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2262:1805:27696:GTGTCAT_GCTAACT\t163\t22\t16345792\t40\t143M\t=\t16346018\t369\tCACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAA\tFFFFFFFFFFF:FFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF:FFFFFFFF:FFFFFFF:,FFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFF::FFFF\tXA:Z:14,-19495406,143M,0;14,+20077911,143M,0;2,-131955560,143M,5;17,+29541346,15M1D128M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGTCAT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2153:24198:24799:GTGTCAT_GCTAACT\t163\t22\t16345792\t40\t143M\t=\t16346018\t369\tCACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGTATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF,FFFFFFFFFFFFFFF,FFF:F\tXA:Z:14,+20077911,143M,1;14,-19495406,143M,1;17,+29541346,15M1D128M,6;\tMC:Z:143M\tMD:Z:67G75\tUI:Z:GTGTCAT_GCTAACT\tNM:i:1\tAS:i:138\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2153:17002:16438:ACAACTC_GCCATTC\t163\t22\t16345793\t40\t143M\t=\t16345948\t298\tACATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFF:F\tXA:Z:14,+20077912,143M,0;14,-19495405,143M,0;2,-131955559,143M,5;17,+29541347,14M1D129M,5;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2125:32479:24298:GTGTCAT_GTGAGAT\t163\t22\t16345795\t40\t143M\t=\t16345933\t281\tATTTGGTAGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF,FF:FFFF,FFFFFFFFFFF:FFFFFFF:FFFFFFFFFFFFFF:FFFFFF:F:F:FFFFFFFFFFF,FFFFFF:F\tXA:Z:14,-19495403,143M,0;14,+20077914,143M,0;2,-131955557,143M,5;17,+29541349,12M1D131M,5;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2123:6840:6856:GTGAGAT_GTCACTC\t99\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF::FFFFFFFFFFF\tXA:Z:14,-19495395,143M,0;14,+20077922,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2152:24442:24314:GTGAGAT_GTCACTC\t99\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495395,143M,0;14,+20077922,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2215:6813:28009:GTGAGAT_GTCACTC\t99\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCAGAATAAATTTATTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFF,FFFFF,FFFFFFFFFF,FFFFFFFFFFFF:FFFFFFFFFFFFFF,FFFFF:F:F:F:FFFFF::,:F:,FF:FFFFF:FFFF,FF,FF,:FFFFFF,F,F\tXA:Z:14,+20077922,143M,2;14,-19495395,143M,2;17,+29541357,4M1D139M,7;2,-131955549,143M,7;\tMC:Z:143M\tMD:Z:32T9C100\tUI:Z:GTGAGAT_GTCACTC\tNM:i:2\tAS:i:133\tXS:i:133",
//                "A01524:289:HFJLTDRX3:1:2243:23863:19053:GTGAGAT_GTCACTC\t99\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF\tXA:Z:14,+20077922,143M,0;14,-19495395,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2243:31358:36010:GTCACTC_GTGAGAT\t163\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFF:FFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF::FFFFFFF,FF\tXA:Z:14,-19495395,143M,0;14,+20077922,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2128:16649:33646:GTCACTC_CTGAGAT\t163\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFF,FFFFFFFFF,FFFFF:FFFFFFFFFFFFF:,FFF,F,FFFF,FFF:FFFFFFFFFFFFFFFFFFF,FFFFFF:FFFFFFFF:FFF,FFFFFFFF:FFF,FF:FFFFFFFFFFFFFFFFFFFFFFFFFFF,:F,\tXA:Z:14,+20077922,143M,0;14,-19495395,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2136:8467:11804:GTGAGAT_GTCACTC\t99\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF::FFFFFFFFFF:FFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFF\tXA:Z:14,-19495395,143M,0;14,+20077922,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2222:6370:1939:TCGTGTC_CAGTACT\t99\t22\t16345803\t0\t24S119M\t=\t16345945\t285\tTATTTTTATACTACCAAATGTGTGGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF,FFFFF:FFFFFFFFFFFFF,FFFFFFF:FFFFF,:FFFFFFFFFF,:FFFFFFFF,FFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFF:FFFFFFF\tSA:Z:14,19495527,+,30M113S,0,0;\tXA:Z:14,-19495419,119M24S,0;14,+20077922,24S119M,0;17,+29541362,28S115M,4;2,-131955573,119M24S,5;\tMC:Z:143M\tMD:Z:119\tUI:Z:TCGTGTC_CAGTACT\tNM:i:0\tAS:i:119\tXS:i:119",
//                "A01524:289:HFJLTDRX3:2:2225:24424:33583:GTGAGAT_GTCACTC\t99\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAATAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFF:FFFFFF,FFFF:FFFFFFFFFFFF:FFFFFFF:FF::FFFFFFFFF,FFFFFFFFFF:FF:FFFFFFFFFFFFFFF:FFFFF:FFFF::FFF:FFFFFFFFFFFFFF:,FFFFFF,FFFFFF,FF\tXA:Z:14,+20077922,143M,1;14,-19495395,143M,1;17,+29541357,4M1D139M,6;2,-131955549,143M,6;\tMC:Z:143M\tMD:Z:126G16\tUI:Z:GTGAGAT_GTCACTC\tNM:i:1\tAS:i:138\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2225:28872:24533:TAGTGTC_CAGTACT\t99\t22\t16345803\t0\t24S119M\t=\t16345945\t285\tTATTTTTATACTACCAAATGTGTGGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF,FFFFFFFFFF:,FFFFF,FF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF:FFFFFFFFFFFFFFFFFFF\tSA:Z:14,19495527,+,30M113S,0,0;\tXA:Z:14,+20077922,24S119M,0;14,-19495419,119M24S,0;17,+29541362,28S115M,4;2,-131955573,119M24S,5;\tMC:Z:143M\tMD:Z:119\tUI:Z:TCGTGTC_CAGTACT\tNM:i:0\tAS:i:119\tXS:i:119",
//                "A01524:289:HFJLTDRX3:2:2232:26169:20290:GTGAGAT_GTCACTC\t99\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFF,,FF\tXA:Z:14,+20077922,143M,0;14,-19495395,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2248:15447:32127:GTCACTC_GTGAGAT\t163\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFF:FFFF:FFFFFFFFFFF\tXA:Z:14,-19495395,143M,0;14,+20077922,143M,0;17,+29541357,4M1D139M,5;2,-131955549,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGAGAT_GTCACTC\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2250:11704:23296:GTCACTC_GTGAGAT\t163\t22\t16345803\t40\t143M\t=\t16345981\t321\tGTATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACACGCCT\tFFF::FFFFF:FFFFFFF:FFFFFFFFF:FF:FFFFFFFFFFFFFFF:FFF:FFFF:F:FFFFFFFFFFFFFFFF:FFFF:FF:FFFF,FFF:F,:FFFFFFFFFFFF,FFFF:FFFFFFFFFFF,FFF:FFFFFFFF,,FFF\tXA:Z:14,-19495395,143M,1;14,+20077922,143M,1;17,+29541357,4M1D139M,6;2,-131955549,143M,6;\tMC:Z:143M\tMD:Z:138A4\tUI:Z:GTGAGAT_GTCACTC\tNM:i:1\tAS:i:138\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2125:30056:10300:CGTGTTT_ACTAGGT\t163\t22\t16345805\t25\t143M\t=\t16346079\t417\tATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATATCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTA\tFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFF:FFFFFFFFFF:F::FFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:F:F:F:FFFFFF:FFF,FFF:F,FF\tXA:Z:14,+20077924,143M,1;14,-19495393,143M,1;17,+29541362,2S141M,5;2,-131955547,143M,6;\tMC:Z:143M\tMD:Z:86C56\tUI:Z:CGTGTTT_ACTAGGT\tNM:i:1\tAS:i:138\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2224:29351:5071:CGTGTTT_ACTAGGT\t163\t22\t16345805\t25\t143M\t=\t16346079\t417\tATAAAAATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATATCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTA\t:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFF:FFF:F,F:FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077924,143M,1;14,-19495393,143M,1;17,+29541362,2S141M,5;2,-131955547,143M,6;\tMC:Z:143M\tMD:Z:86C56\tUI:Z:CGTGTTT_ACTAGGT\tNM:i:1\tAS:i:138\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2266:12988:23829:GTGTCAT_CGATATG\t99\t22\t16345811\t40\t143M\t=\t16346048\t380\tATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFF:FFFFFFFFFFF\tXA:Z:14,+20077930,143M,0;14,-19495387,143M,0;2,-131955541,143M,4;17,+29541366,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGTCAT_CGATATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2155:24930:23594:GTGTCAT_CGATATG\t99\t22\t16345811\t40\t143M\t=\t16346048\t380\tATAGCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF\tXA:Z:14,+20077930,143M,0;14,-19495387,143M,0;2,-131955541,143M,4;17,+29541366,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:GTGTCAT_CGATATG\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2111:24939:2848:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFF:FFF,:FFF:F,:,:FFFFFFF,:FFFFFFFFFFFF,FF\tXA:Z:14,-19495383,143M,0;14,+20077934,143M,0;17,+29541370,143M,4;2,-131955537,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2128:5159:25927:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF,FF,:FFFFFFF:FFF::F\tXA:Z:14,+20077934,143M,0;14,-19495383,143M,0;17,+29541370,143M,4;2,-131955537,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2142:28637:34491:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFF:FFFFFFFF\tXA:Z:14,-19495383,143M,0;14,+20077934,143M,0;17,+29541370,143M,4;2,-131955537,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2158:17155:11256:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFF:FFFF:,FFFFFFFF,FFFFFFFF::FFFFFFFFFF:FFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFF:FFFFF,:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:,FFFFFFFFFFF,FFF::F\tXA:Z:14,+20077934,143M,0;14,-19495383,143M,0;17,+29541370,143M,4;2,-131955537,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2170:25690:33990:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495383,143M,0;14,+20077934,143M,0;2,-131955537,143M,4;17,+29541370,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2272:5141:8516:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGACTTAAATTTAAAGA\tFF:FFFFFFFF:FFFF:::FFFF:FFFFFF:FFFF::F:FF:FFFFFF::FFFFFF:FFFFFFFFFFFF:FFFFFFF,FFF:FFFFF:FFFFFFFFFFFF::FFFFFFFFFFFF:FF:F,:FFFFF:,,FFFF:FFFFF,:FF\tXA:Z:14,+20077934,143M,1;14,-19495383,143M,1;17,+29541370,143M,5;2,-131955537,143M,5;\tMC:Z:143M\tMD:Z:128C14\tUI:Z:CAGACGT_ATCTCCT\tNM:i:1\tAS:i:138\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2273:19687:4523:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FF,FFFFFFFFFFFFFFFFFFFFFFFFF:FFFF,FFFFF:F:F\tXA:Z:14,+20077934,143M,0;14,-19495383,143M,0;2,-131955537,143M,4;17,+29541370,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2217:23918:8970:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077934,143M,0;14,-19495383,143M,0;17,+29541370,143M,4;2,-131955537,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2245:30680:13385:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFFFFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495383,143M,0;14,+20077934,143M,0;2,-131955537,143M,4;17,+29541370,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2254:28655:16736:CAGACGT_ATCTCCT\t163\t22\t16345815\t40\t143M\t=\t16345961\t289\tCTTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGA\tFFFFFFFFFFFFFFFFFFF:FFF::FFFFFFFFFFF:F:FFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFF,FFFF:,FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF,FF:FFFFFFF:FFF:FFFF:FF\tXA:Z:14,-19495383,143M,0;14,+20077934,143M,0;17,+29541370,143M,4;2,-131955537,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:CAGACGT_ATCTCCT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2258:1922:9142:TGTCGTG_TCGTGTT\t163\t22\t16345816\t40\t13S130M\t=\t16346048\t375\tTTGAAAAAAGCTATTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFFFFF:FF:FFFFFFFFFFF,FFFF:F:FFFFFFFFFFFFFFFFFFF:FFFFFFFFF,:,FFFFFFFFFF:F,FFF\tXA:Z:14,-19495395,130M13S,0;14,+20077935,13S130M,0;17,+29541371,13S130M,4;2,-131955549,130M13S,4;\tMC:Z:143M\tMD:Z:130\tUI:Z:TGTCGTG_TCGTGTT\tNM:i:0\tAS:i:130\tXS:i:130",
//                "A01524:289:HFJLTDRX3:1:2264:11650:34100:TGTCGTG_TCGTGTT\t163\t22\t16345816\t40\t13S130M\t=\t16346048\t375\tTTGAAAAAAGCTATTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFF::FFFFFFFFFFFFFFFF:FFFFFFFFFF,FFFFF,FFFFFFF:FFFFF:FFF:FFFFFFFFFFFFFF,FFFF:FFFFFFFF,F:FFFFFF:FFFF:FFFFFFF:FFFFFFFFFFFF,F:FFF:F:FFFFFFF,FF\tXA:Z:14,+20077935,13S130M,0;14,-19495395,130M13S,0;2,-131955549,130M13S,4;17,+29541371,13S130M,4;\tMC:Z:143M\tMD:Z:130\tUI:Z:TGTCGTG_TCGTGTT\tNM:i:0\tAS:i:130\tXS:i:130",
//                "A01524:289:HFJLTDRX3:2:2210:6126:28729:TGTCGTG_TCGTGTT\t163\t22\t16345816\t40\t13S130M\t=\t16346048\t375\tTTGAAAAAAGCTATTTTTTCAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:F:FFF:FFFF:FF,FFF\tXA:Z:14,-19495395,130M13S,0;14,+20077935,13S130M,0;2,-131955549,130M13S,4;17,+29541371,13S130M,4;\tMC:Z:143M\tMD:Z:130\tUI:Z:TGTCGTG_TCGTGTT\tNM:i:0\tAS:i:130\tXS:i:130",
//                "A01524:289:HFJLTDRX3:1:2233:20473:18505:CATGATA_ACAACTC\t163\t22\t16345823\t40\t143M\t=\t16345955\t275\tAAAGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFF:FFF:F:FFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077942,143M,0;14,-19495375,143M,0;17,+29541378,143M,4;2,-131955529,143M,4;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2143:6370:34976:AGCATCT_GCTAACT\t99\t22\t16345826\t40\t143M\t=\t16345930\t247\tGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF:F:FFFF:FF:FFF\tXA:Z:14,+20077945,143M,0;14,-19495372,143M,0;17,+29541389,8S135M,1;2,-131955526,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2258:22444:10614:AGCATCT_GCTAACT\t99\t22\t16345826\t40\t143M\t=\t16345930\t247\tGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF\tXA:Z:14,-19495372,143M,0;14,+20077945,143M,0;17,+29541389,8S135M,1;2,-131955526,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:1:2273:32687:5995:AGCATCT_GCTAACT\t99\t22\t16345826\t40\t143M\t=\t16345930\t247\tGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAG\tFFFFFFFFFFF,FFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF:FFF:FF\tXA:Z:14,-19495372,143M,0;14,+20077945,143M,0;17,+29541389,8S135M,1;2,-131955526,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2121:21269:31845:AGCATCT_GCTAACT\t99\t22\t16345826\t40\t143M\t=\t16345930\t247\tGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFF\tXA:Z:14,+20077945,143M,0;14,-19495372,143M,0;17,+29541389,8S135M,1;2,-131955526,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2212:6488:27853:AGCATCT_GCTAACT\t99\t22\t16345826\t40\t143M\t=\t16345930\t247\tGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAG\tFF,FF:FFFFFFFFFFFFF:FFFFFFF::FFF:FFFFFFF,F::FFFFFFFFFFFFF:,FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:,,FFFFFFFFFFFFFFFFF,:FFFF\tXA:Z:14,-19495372,143M,0;14,+20077945,143M,0;17,+29541389,8S135M,1;2,-131955526,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2212:7898:28166:AGCATCT_GCTAACT\t99\t22\t16345826\t40\t143M\t=\t16345930\t247\tGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFFFF,F,:FFFFFFFFFF\tXA:Z:14,+20077945,143M,0;14,-19495372,143M,0;17,+29541389,8S135M,1;2,-131955526,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2222:19985:28557:AGCATCT_GCTAACT\t99\t22\t16345826\t40\t143M\t=\t16345930\t247\tGTTTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,:FFFFF\tXA:Z:14,-19495372,143M,0;14,+20077945,143M,0;17,+29541389,8S135M,1;2,-131955526,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_GCTAACT\tNM:i:0\tAS:i:143\tXS:i:143",
//                "A01524:289:HFJLTDRX3:2:2139:9471:19179:AATGCCT_TCGTGTC\t99\t22\t16345829\t40\t143M\t=\t16346026\t340\tTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCT\tFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFF:FFFFFF,FFFFFFFFFFFFF::FFF,FFFFF:FFFF\tXA:Z:14,-19495369,143M,1;14,+20077948,143M,1;17,+29541389,5S135M3S,1;2,-131955526,3S140M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_TCGTGTC\tNM:i:0\tAS:i:143\tXS:i:142",
//                "A01524:289:HFJLTDRX3:2:2140:16957:22780:AATGCCT_TCGTGTC\t99\t22\t16345829\t40\t143M\t=\t16346026\t340\tTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCT\tFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,:FF:FFFFFFFFF,FFFFFFFFFFFFFFF\tXA:Z:14,-19495369,143M,1;14,+20077948,143M,1;17,+29541389,5S135M3S,1;2,-131955526,3S140M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_TCGTGTC\tNM:i:0\tAS:i:143\tXS:i:142",
//                "A01524:289:HFJLTDRX3:2:2166:7419:29684:AATGCCT_TCGTGTC\t99\t22\t16345829\t40\t143M\t=\t16346026\t340\tTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFF:FFFFF:FFF\tXA:Z:14,-19495369,143M,1;14,+20077948,143M,1;17,+29541389,5S135M3S,1;2,-131955526,3S140M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_TCGTGTC\tNM:i:0\tAS:i:143\tXS:i:142",
//                "A01524:289:HFJLTDRX3:2:2207:13250:33771:AATGCCT_TCGTGTC\t99\t22\t16345829\t40\t143M\t=\t16346026\t340\tTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFF,F,F:FFF:,FFF::\tXA:Z:14,-19495369,143M,1;14,+20077948,143M,1;17,+29541389,5S135M3S,1;2,-131955526,3S140M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_TCGTGTC\tNM:i:0\tAS:i:143\tXS:i:142",
//                "A01524:289:HFJLTDRX3:2:2210:7139:9940:AATGCCT_TCGTGTC\t99\t22\t16345829\t40\t143M\t=\t16346026\t340\tTTCTCATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:F:FFFFFFF:FFFFFFFF:FFF:,FFFFFFFFFF:FFFFFFFFFFFFFF:,FFFFF:FFF:FF,FFFFFFFFFFFFF,FFFF::FFF,F:F:F\tXA:Z:14,+20077948,143M,1;14,-19495369,143M,1;17,+29541389,5S135M3S,1;2,-131955526,3S140M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_TCGTGTC\tNM:i:0\tAS:i:143\tXS:i:142",
//                "A01524:289:HFJLTDRX3:1:2158:10041:5321:GCCATTG_TTGGCTC\t163\t22\t16345834\t46\t143M\t=\t16346016\t325\tATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAG\tFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:F,FFF:F:FFFFFFFFFFFFF\tXA:Z:14,-19495364,143M,1;14,+20077953,143M,1;17,+29541389,143M,3;2,-131955518,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GCCATTG_TTGGCTC\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2169:13196:35368:GCCATTG_TTGGCTC\t163\t22\t16345834\t46\t143M\t=\t16346016\t325\tATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFF:F:FFFFFFFFF:FFFFFFFFF:FFFFFFF:,FFFFFFFF::FFFFF\tXA:Z:14,-19495364,143M,1;14,+20077953,143M,1;17,+29541389,143M,3;2,-131955518,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GCCATTG_TTGGCTC\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2219:3794:34929:GCCATTG_TTGGCTC\t163\t22\t16345834\t46\t143M\t=\t16346016\t325\tATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAG\tFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF:FFFFFFFFF:FFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF\tXA:Z:14,+20077953,143M,1;14,-19495364,143M,1;17,+29541389,143M,3;2,-131955518,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GCCATTG_TTGGCTC\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2265:13449:1454:GCCATTG_TTGGCTC\t163\t22\t16345834\t46\t143M\t=\t16346016\t325\tATAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAG\tFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF::F:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFF:F:F:F:,FFF,,FF\tXA:Z:14,+20077953,143M,1;14,-19495364,143M,1;17,+29541389,143M,3;2,-131955518,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:GCCATTG_TTGGCTC\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2124:8684:23109:CAGTACT_CTTGGAT\t163\t22\t16345836\t46\t143M\t=\t16345911\t218\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF::FFF:FFFF\tXA:Z:14,-19495362,143M,1;14,+20077955,143M,1;17,+29541391,143M,3;2,-131955516,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CTTGGAT_CAGTACT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2144:15320:18818:AAGGATA_TGTCGTT\t163\t22\t16345836\t46\t143M\t=\t16345919\t226\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF:FFFFFFF,F:FF\tXA:Z:14,+20077955,143M,1;14,-19495362,143M,1;17,+29541391,143M,3;2,-131955516,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:AAGGATA_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2144:16468:18333:AAGGATA_TGTCGTT\t163\t22\t16345836\t46\t143M\t=\t16345919\t226\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGAC\tFFFFFFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFF:F::FFFFFFFFFFFFFFFFFFFFFFF:FFF:,FFFFFF:FFFFFFFFFFFFF,FFFF:,FFF,FF::F:FFFFF:\tXA:Z:14,-19495362,143M,1;14,+20077955,143M,1;17,+29541391,143M,3;2,-131955516,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:AAGGATA_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2219:24777:32189:AAGGATA_TGTCGTT\t163\t22\t16345836\t46\t143M\t=\t16345919\t226\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFF\tXA:Z:14,+20077955,143M,1;14,-19495362,143M,1;17,+29541391,143M,3;2,-131955516,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:AAGGATA_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2227:28483:10144:AAGGATA_TGTCGTT\t163\t22\t16345836\t46\t143M\t=\t16345919\t226\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTAAAGGCTTGGAGAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF,FFFFFFF,FFFFFFFFFFFFFFFF:FF,FFFFFFFFFFF::FFF:FFFFFFFFFFFF,FFFFF:F:F,F,FFFFF,:,FF\tXA:Z:14,-19495362,143M,2;14,+20077955,143M,2;17,+29541391,143M,4;2,-131955516,143M,6;\tMC:Z:143M\tMD:Z:130C12\tUI:Z:AAGGATA_TGTCGTT\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:1:2227:29776:8813:AAGGATA_TGTCGTT\t163\t22\t16345836\t46\t143M\t=\t16345919\t226\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFF,FFFFF:FFF:FFFF:\tXA:Z:14,+20077955,143M,1;14,-19495362,143M,1;17,+29541391,143M,3;2,-131955516,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:AAGGATA_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2101:5502:17002:AAGGATA_TGTCGTT\t163\t22\t16345836\t46\t143M\t=\t16345919\t226\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATTTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGAC\tFFFFFFFF:FF,,FFFF:F,F:,:FFFFFFF,:F,FFFFFFF:FFF:,:FFFF:F,:F:,,FF:F:FFFFFFF::FFFF,FFFF:FFFF,FF:,FFFFFF,FFFFF,,F,,FF:,,FFFF,:FFFFFF,FFF,,FF,,::FFF\tXA:Z:14,+20077955,143M,2;14,-19495362,143M,2;17,+29541391,143M,4;2,-131955516,143M,6;\tMC:Z:143M\tMD:Z:59G83\tUI:Z:AAGGATA_TGTCGTT\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:2:2238:8386:29857:CTTGGAT_CAGTACT\t99\t22\t16345836\t46\t143M\t=\t16345911\t218\tAATAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495362,143M,1;14,+20077955,143M,1;17,+29541391,143M,3;2,-131955516,143M,5;\tMC:Z:143M\tMD:Z:143\tUI:Z:CTTGGAT_CAGTACT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2142:22851:20964:GTCGTTG_CGTGTTA\t163\t22\t16345838\t46\t143M\t=\t16346057\t362\tTAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACATACAGGCTTGGAGACAA\tFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FF:FFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,:,FFFFFFF:F:FF,FFF:F,FF:F::FF::F:FFFFF\tXA:Z:14,-19495360,143M,2;14,+20077957,143M,2;17,+29541393,143M,4;2,-131955514,143M,7;\tMC:Z:143M\tMD:Z:125C17\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:1:2161:2193:35446:CGGTTGT_TCGTGTA\t163\t22\t16345838\t46\t143M\t=\t16345970\t275\tTAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFF:FFFF:FFFFFFFF::FFFFFFFFFFF:FFFFF:FFFF\tXA:Z:14,-19495360,143M,1;14,+20077957,143M,1;17,+29541393,143M,3;2,-131955514,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCGTGTA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2162:3125:5337:CGGTTGT_TCGTGTA\t163\t22\t16345838\t46\t143M\t=\t16345970\t275\tTAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FFF:FFFFFFFFFFFF:FFFF:FFFFFFF:FFFFFFFFFF\tXA:Z:14,-19495360,143M,1;14,+20077957,143M,1;17,+29541393,143M,3;2,-131955514,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCGTGTA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2221:4399:29121:CGGTTGT_TCGTGTA\t163\t22\t16345838\t46\t143M\t=\t16345970\t275\tTAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAA\tFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFF:F:FFFFFFF:FFFFFFFFF:FFF,:FFFF,F:FFFF,FFF::FFF,FF\tXA:Z:14,+20077957,143M,1;14,-19495360,143M,1;17,+29541393,143M,3;2,-131955514,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCGTGTA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2222:13268:25316:CGGTTGT_TCGTGTA\t163\t22\t16345838\t46\t143M\t=\t16345970\t275\tTAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAA\tFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FF:FF,FFFFFFFFFFF\tXA:Z:14,+20077957,143M,1;14,-19495360,143M,1;17,+29541393,143M,3;2,-131955514,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCGTGTA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2177:28049:30060:CGGTTGT_TCGTGTA\t163\t22\t16345838\t46\t143M\t=\t16345970\t275\tTAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAA\tFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFF:FFFFFFFFFFFF,FFFF:FFFFFF\tXA:Z:14,-19495360,143M,1;14,+20077957,143M,1;17,+29541393,143M,3;2,-131955514,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCGTGTA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2234:14904:8328:CGGTTGT_TCGTGTA\t163\t22\t16345838\t46\t143M\t=\t16345970\t275\tTAAATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFF:FF\tXA:Z:14,-19495360,143M,1;14,+20077957,143M,1;17,+29541393,143M,3;2,-131955514,143M,6;\tMC:Z:143M\tMD:Z:143\tUI:Z:CGGTTGT_TCGTGTA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2139:21911:3709:AACACAT_CATGATA\t163\t22\t16345841\t46\t143M\t=\t16345940\t242\tATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAATAAAAACCTACAGGCTTGGAGACAAGAA\tFFFF:FFFFFFFFFFF:FFFFFFFFFF,FFFFFFFFFFFFFFFFFFFF:FF:F:FFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFF:F:FFFF,FFF,FF,,:FFFFFFFFFF:,FFFFFFF:F,F,F:FF,FF,:,FF:FF\tXA:Z:14,-19495357,143M,2;14,+20077960,143M,2;17,+29541396,143M,4;\tMC:Z:143M\tMD:Z:115G27\tUI:Z:AACACAT_CATGATA\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:2:2248:10176:17707:AACACAT_CATGATA\t163\t22\t16345841\t46\t143M\t=\t16345940\t242\tATTTCTTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495357,143M,1;14,+20077960,143M,1;17,+29541396,143M,3;\tMC:Z:143M\tMD:Z:143\tUI:Z:AACACAT_CATGATA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2110:21748:30076:GCCATTG_AATGCCT\t163\t22\t16345847\t46\t143M\t=\t16345951\t247\tTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:F,FFFF\tXA:Z:14,-19495351,143M,1;14,+20077966,143M,1;17,+29541402,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_GACATTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2236:31982:3051:GCCATTG_AATGCCT\t163\t22\t16345847\t45\t143M\t=\t16345951\t247\tTTTAAAAGGTTGGATAGCTATTATCCTCAGTCTTATGTCTGATACCATGATTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACATAAATGCTTGGAGACAAGAAGCTGTA\t:FFFFFF,,FFF,FFF:F,FFF:F:FF,FF:F,FF,FFFF,F,FFF:FF,F,FFF:FFFFFF,:FFFF::F,FFFFF,FF:FFFFFFFFFFFFFF,FFFFFFFF,FFFFFFFFFF:,:F:F,F,FFFFFFFFF:,FF:,F,FF\tXA:Z:14,-19495351,143M,6;14,+20077966,143M,6;17,+29541402,143M,9;\tMC:Z:143M\tMD:Z:27G21T66C2C1G21\tUI:Z:AATGCCT_GACATTG\tNM:i:5\tAS:i:118\tXS:i:113",
//                "A01524:289:HFJLTDRX3:2:2126:21920:5040:AATGCCT_GACATTG\t99\t22\t16345847\t46\t143M\t=\t16345951\t247\tTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTA\tFF:FFF,FF::F,FFF:,FF:FFFFF::F:F::,F:,FFF,F:FF:FF,F,FFFFF,F:,FF,:FF,F:FF:,:,FF,,FFFFFFFFF,:FF,FF,F:F,FFF:FFFFF:FFFF,F:FF,FF,FF,FFF,FFFF,FFFF:FFF\tXA:Z:14,-19495351,143M,1;14,+20077966,143M,1;17,+29541402,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_GACATTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2218:27380:19476:GCCATTG_AATGCCT\t163\t22\t16345847\t46\t143M\t=\t16345951\t247\tTTTAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTA\tFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFF,FFFFFF:FF::FFF:F::FF:FFFFF\tXA:Z:14,-19495351,143M,1;14,+20077966,143M,1;17,+29541402,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AATGCCT_GACATTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2253:14733:20181:TTGGCTT_TGTCGTT\t99\t22\t16345850\t46\t143M\t=\t16345892\t185\tAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:F:FFFFFF:F:FFFFFFF\tXA:Z:14,+20077969,143M,1;14,-19495348,143M,1;17,+29541405,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:TTGGCTT_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2148:4788:34491:TTGGCTT_TGTCGTT\t99\t22\t16345850\t46\t143M\t=\t16345892\t185\tAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495348,143M,1;14,+20077969,143M,1;17,+29541405,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:TTGGCTT_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2266:19153:33880:TTGGCTT_TGTCGTT\t99\t22\t16345850\t46\t143M\t=\t16345892\t185\tAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGT\tFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFF:FFFFF,,FFFFFFFFFFF\tXA:Z:14,+20077969,143M,1;14,-19495348,143M,1;17,+29541405,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:TTGGCTT_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2107:6470:28197:ATCTCCT_TGTCGTC\t99\t22\t16345854\t46\t143M\t=\t16345893\t182\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495344,143M,1;14,+20077973,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ATCTCCT_TGTCGTC\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2116:28971:5697:ATCTCCT_CACTGTG\t99\t22\t16345854\t46\t143M\t=\t16346030\t319\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ATCTCCT_CACTGTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2217:1642:34397:AGCATCT_AAAGATA\t163\t22\t16345854\t40\t143M\t=\t16345975\t264\tGGTTGGATAGCTATTATACTGAGTCTTATGTCTGATACCATGTTTTTTTTTTGTTTTTAGAGTCTTATATTAAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGATTCTATACAAGAATCTGTAAGTATCT\t,,F:,,F,FFF,FFFF,,F,,F,F,,FFFFFF:FF,FF,F,FFFFFF,:FF::FFF:FFFFFF,,FF:FFF,FFFFFFFF:,FF:FF::F,,:FFF:F:F:F,,FF:F,::F,F,,F:F,:FFF:FFFFF,F:,,FF:F,:,F\tXA:Z:14,+20077973,143M,9;14,-19495344,143M,9;17,+29541409,115M28S,4;2,-131955526,28S115M,6;\tMC:Z:143M\tMD:Z:17C29G23T44C2G0G1G7G12\tUI:Z:AGCATCT_CATGATA\tNM:i:8\tAS:i:103\tXS:i:101",
//                "A01524:289:HFJLTDRX3:1:2241:17897:33520:ATCTCCT_CACTGTG\t99\t22\t16345854\t46\t143M\t=\t16346030\t319\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ATCTCCT_CACTGTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2262:7075:22263:AGCATCT_CATGATA\t163\t22\t16345854\t46\t143M\t=\t16345975\t264\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF:FFFFFFFFF:FFFFF:F:FFFFFF,,:FFFFF:FFFF:F:F:FFFFFFFF::F,FFFFFFFFFFFF\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_CATGATA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2122:31485:32221:AGCATCT_CATGATA\t163\t22\t16345854\t46\t143M\t=\t16345975\t264\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495344,143M,1;14,+20077973,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_CATGATA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2223:16957:6621:GTGTCAT_AGCATCT\t163\t22\t16345854\t46\t143M\t=\t16345892\t181\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF:FFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495344,143M,1;14,+20077973,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2234:30364:31469:AGCATCT_CATGATA\t163\t22\t16345854\t46\t143M\t=\t16345975\t264\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFF,FFFFFFF:FFFFFFFF,FFFFFF\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_CATGATA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2246:20916:7404:AGCATCT_CATGATA\t163\t22\t16345854\t46\t143M\t=\t16345975\t264\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFF:FFFFFFFFF:F\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_CATGATA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2264:27127:21699:TGTCGTC_ATCTCCT\t163\t22\t16345854\t46\t143M\t=\t16345893\t182\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:,FFFFFFFFF:F::FFFFFFFF:FFFFF:FFFFFFFFFFFF\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ATCTCCT_TGTCGTC\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2271:31304:36260:AGCATCT_CATGATA\t163\t22\t16345854\t46\t143M\t=\t16345975\t264\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\t:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFF,FFFF:FFFFFFFFFF:FFFFFFFFFFFF\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:AGCATCT_CATGATA\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2277:27290:19132:ATCTCCT_CACTGTT\t99\t22\t16345854\t46\t143M\t=\t16346030\t319\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFF,:FF,:FF,FFFFFFFFFFFF:FFFF,FFF:F,FF,FF:FFFF,F:F::F,,FFFFFF:FFFFFF,,F,F:F,FFFFFFFFF,FFFFFFFFF:FF:,:FFF:FFFF,F,FFFF,FF,FFF::F:,FFFFFFF,FFF,\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ATCTCCT_CACTGTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2235:5502:16501:ACTAGGT_TCGTGTG\t163\t22\t16345855\t46\t143M\t=\t16345890\t178\tGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCTT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFF:F:FFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF\tXA:Z:14,+20077974,143M,1;14,-19495343,143M,1;17,+29541410,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ACTAGGT_TCGTGTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2117:30074:23422:ACTAGGT_TCGTGTG\t163\t22\t16345855\t46\t143M\t=\t16345890\t178\tGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCTT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495343,143M,1;14,+20077974,143M,1;17,+29541410,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ACTAGGT_TCGTGTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2158:29080:15186:TTGGCTC_GTCGTTG\t163\t22\t16345855\t46\t143M\t=\t16345899\t187\tGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCTT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFF:FFFFF,FFFFFFFFFFFF:FFF:FFFF:FFFFFFF::FFFFF:FFFF\tXA:Z:14,+20077974,143M,1;14,-19495343,143M,1;17,+29541410,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:TTGGCTC_GTCGTTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2237:6714:36323:TTGGCTC_GTCGTTG\t163\t22\t16345855\t46\t143M\t=\t16345899\t187\tGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCTT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF:FFFFFF,F:FFFFFFFFFFF,FFFFFFFFFFF:FFFF:FFFFFFF,FF,,FFFFF,FFFFFF\tXA:Z:14,+20077974,143M,1;14,-19495343,143M,1;17,+29541410,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:TTGGCTC_GTCGTTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:2:2261:8323:27868:ACTAGGT_TCGTGTG\t163\t22\t16345855\t46\t143M\t=\t16345890\t178\tGTTGGTTAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCTT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495343,143M,2;14,+20077974,143M,2;\tMC:Z:143M\tMD:Z:5A137\tUI:Z:ACTAGGT_TCGTGTG\tNM:i:1\tAS:i:138\tXS:i:133",
//                "A01524:289:HFJLTDRX3:2:2270:4571:2895:TTGGCTC_GTCGTTG\t163\t22\t16345855\t46\t143M\t=\t16345899\t187\tGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCTT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFF:FFFF:FFFF:FFFFFFFF:F:FFFFFFF,FFFFFFFFFFFF\tXA:Z:14,+20077974,143M,1;14,-19495343,143M,1;17,+29541410,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:TTGGCTC_GTCGTTG\tNM:i:0\tAS:i:143\tXS:i:138",
//                "A01524:289:HFJLTDRX3:1:2210:31955:26741:AACACAT_TACGATA\t83\t22\t16345893\t48\t42S101M\t=\t16345894\t-100\tGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTACGATATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-20078012,42S101M,1;14,+19495347,101M42S,1;\tMC:Z:100M43S\tMD:Z:101\tUI:Z:AACACAT_TACGATA\tNM:i:0\tAS:i:101\tXS:i:96",
//                "A01524:289:HFJLTDRX3:1:2210:32090:26412:AACACAT_TACGATA\t83\t22\t16345893\t48\t42S101M\t=\t16345894\t-100\tGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTACGATATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,+19495347,101M42S,1;14,-20078012,42S101M,1;\tMC:Z:100M43S\tMD:Z:101\tUI:Z:AACACAT_TACGATA\tNM:i:0\tAS:i:101\tXS:i:96",
//                "A01524:289:HFJLTDRX3:2:2269:10963:11960:AACACAT_TACGATA\t83\t22\t16345893\t48\t42S101M\t=\t16345894\t-100\tGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTACGATATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-20078012,42S101M,1;14,+19495347,101M42S,1;\tMC:Z:100M43S\tMD:Z:101\tUI:Z:AACACAT_TACGATA\tNM:i:0\tAS:i:101\tXS:i:96"
//                );
//
//        boolean finished = false;
//        while(!finished)
//        {
//            finished = true;
//            for(int i = 0; i < samStrings.size(); i++)
//            {
//                final int idx = i;
//                List<String> subset = IntStream.range(0, samStrings.size()).filter(j -> j != idx).mapToObj(j -> samStrings.get(j)).collect(Collectors.toList());
//                if(!testSubsetOfReads(subset))
//                {
//                    samStrings.clear();
//                    samStrings.addAll(subset);
//                    finished = false;
//                    break;
//                }
//            }
//        }
//
//        for(String samString0 : samStrings)
//        {
//            String samString = samString0;
//            samString = samString.replaceAll("[\\t]", quoteReplacement("\\t"));
//            samString = format("\"%s\\n\",", samString);
//            System.out.println(samString);
//        }
//
//        // TODO:
////        assertEquals(2, writer.nonConsensusWriteCount());
//        // TODO:
////        assertEquals(1, writer.consensusWriteCount());
//
//        // TODO: remove
//        fail();
//    }
//
//    // TODO: remove this
//    private static boolean testSubsetOfReads(final List<String> samStrings)
//    {
//        String chromosome = "22";
//        int chromosomeLength = 18_000_000;
//
//        MockRefGenome refGenome = new MockRefGenome(true);
//        refGenome.RefGenomeMap.put(chromosome, "A".repeat(chromosomeLength));
//        refGenome.ChromosomeLengths.put(chromosome, chromosomeLength);
//
//        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
//        TestBamWriter writer = new TestBamWriter(config);
//        PartitionReader partitionReader = createPartitionRead(config, writer);
//
//        partitionReader.setupRegion(new ChrBaseRegion(chromosome, 1, chromosomeLength));
//
//        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
//        SAMFileHeader samFileHeader = recordBuilder.getHeader();
//        samFileHeader.setSequenceDictionary(SAM_DICTIONARY_V37);
//        SAMLineParser samLineParser = new SAMLineParser(samFileHeader);
//        for(String samString : samStrings)
//        {
//            partitionReader.processRead(samLineParser.parseLine(samString));
//        }
//
//        try
//        {
//            partitionReader.postProcessRegion();
//        }
//        catch(RuntimeException e)
//        {
//            return false;
//        }
//
//        return true;
//    }

    // TODO: uncomment
//    @Test
//    public void testIlluminaPolyGDuplexUmiGroupCollapse()
//    {
//        MockRefGenome refGenome = new MockRefGenome(true);
//        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
//        refGenome.ChromosomeLengths.put(CHR_1, 1_000);
//
//        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
//        TestBamWriter writer = new TestBamWriter(config);
//        PartitionReader partitionReader = createPartitionRead(config, writer);
//
//        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));
//
//        String umidIdPart1 = "TCCTATG";
//        String umidId1Part2 = "CGGGGGG";
//        String umidId2Part2 = "GGGGGGG";
//        String umiId1 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
//        String umiId2 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;
//
//        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
//                nextReadId(umiId1), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
//        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
//                nextReadId(umiId2), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
//        read2.setMateUnmappedFlag(true);
//
//        partitionReader.processRead(read1);
//        partitionReader.processRead(read2);
//        partitionReader.postProcessRegion();
//
//        assertEquals(2, writer.nonConsensusWriteCount());
//        assertEquals(1, writer.consensusWriteCount());
//    }
//
//    @Test
//    public void testIlluminaPolyGDuplexUmiGroupNoJitterCollapse()
//    {
//        MockRefGenome refGenome = new MockRefGenome(true);
//        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
//        refGenome.ChromosomeLengths.put(CHR_1, 1_000);
//
//        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
//        TestBamWriter writer = new TestBamWriter(config);
//        PartitionReader partitionReader = createPartitionRead(config, writer);
//
//        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));
//
//        String umidIdPart1 = "TCCTATG";
//        String umidId1Part2 = "CGGGGGG";
//        String umidId2Part2 = "GGGGGGG";
//        String umiId1 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
//        String umiId2 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;
//
//        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
//                nextReadId(umiId1), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
//        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
//                nextReadId(umiId2), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
//        read2.setMateUnmappedFlag(true);
//
//        partitionReader.processRead(read1);
//        partitionReader.processRead(read2);
//        partitionReader.postProcessRegion();
//
//        assertEquals(2, writer.nonConsensusWriteCount());
//        assertEquals(0, writer.consensusWriteCount());
//    }
//
//    @Test
//    public void testIlluminaPolyGDuplexUmiGroupNoCollapseWithoutPolyGTail()
//    {
//        MockRefGenome refGenome = new MockRefGenome(true);
//        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
//        refGenome.ChromosomeLengths.put(CHR_1, 1_000);
//
//        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
//        TestBamWriter writer = new TestBamWriter(config);
//        PartitionReader partitionReader = createPartitionRead(config, writer);
//
//        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));
//
//        String umidIdPart1 = "TCCTATG";
//        String umidId1Part2 = "CGGGGGA";
//        String umidId2Part2 = "GGGGGGA";
//        String umiId1 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
//        String umiId2 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;
//
//        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
//                nextReadId(umiId1), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
//        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
//                nextReadId(umiId2), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
//        read2.setMateUnmappedFlag(true);
//
//        partitionReader.processRead(read1);
//        partitionReader.processRead(read2);
//        partitionReader.postProcessRegion();
//
//        assertEquals(2, writer.nonConsensusWriteCount());
//        assertEquals(0, writer.consensusWriteCount());
//    }

//    @Test
//    public void testIlluminaPolyGDuplexUmiGroupCollapseWithJitterCollapsedGroup()
//    {
//        MockRefGenome refGenome = new MockRefGenome(true);
//        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
//        refGenome.ChromosomeLengths.put(CHR_1, 1_000);
//
//        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
//
//        String umidId1Part1 = "TCCTATG";
//        String umidId4Part1 = "TCCTATT";
//        String umidId1Part2 = "GGGGGGG";
//        String umidId2Part2 = "CGGGGGG";
//        String umidId4Part2 = "GGGGGGT";
//
//        String umiId1 = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
//        String umiId2 = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;
//        String umiId3 = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidId4Part2;
//        String umiId4 = umidId4Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
//        String umiId5 = umiId2;
//        String umiId6 = umiId2;
//
//        String readName1 = nextReadId(umiId1);
//        String readName2 = nextReadId(umiId2);
//        String readName3 = nextReadId(umiId3);
//        String readName4 = nextReadId(umiId4);
//        String readName5 = nextReadId(umiId5);
//        String readName6 = nextReadId(umiId6);
//
//        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
//                readName1, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
//        FragmentCoords coords1 = FragmentCoords.fromRead(read1, true);
//        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
//                readName2, CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
//        read2.setMateUnmappedFlag(true);
//        SAMRecord read3 = SamRecordTestUtils.createSamRecord(
//                readName3, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
//        read3.setMateUnmappedFlag(true);
//        SAMRecord read4 = SamRecordTestUtils.createSamRecord(
//                readName4, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
//        read4.setMateUnmappedFlag(true);
//        SAMRecord read5 = SamRecordTestUtils.createSamRecord(
//                readName5, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
//        read5.setMateUnmappedFlag(true);
//        SAMRecord read6 = SamRecordTestUtils.createSamRecord(
//                readName6, CHR_1, 100, TEST_READ_BASES, NO_CIGAR, CHR_1, 100, false, false, null, true, TEST_READ_CIGAR);
//        read6.setReadUnmappedFlag(true);
//
//        List<DuplicateGroup> umiGroups = Lists.newArrayList(
//                new DuplicateGroup(Lists.newArrayList(read2, read3, read4, read5), FragmentCoords.fromRead(read5, true)));
//        List<ReadInfo> singleFragments = Lists.newArrayList(
//                new ReadInfo(read1, coords1), new ReadInfo(read6, FragmentCoords.fromRead(read6, true)));
//        collapsePolyGDuplexUmis(ILLUMINA, config.UMIs, umiGroups, singleFragments);
//
//        assertEquals(1, umiGroups.size());
//        assertEquals(1, singleFragments.size());
//
//        DuplicateGroup umiGroup = umiGroups.get(0);
//        assertEquals(coords1, umiGroup.fragmentCoordinates());
//
//        List<String> expectedCollapsedReadNames = Lists.newArrayList(readName1, readName2, readName3, readName4, readName5, readName6);
//        Collections.sort(expectedCollapsedReadNames);
//        List<String> actualCollapsedReadNames = umiGroup.reads().stream()
//                .map(SAMRecord::getReadName)
//                .collect(Collectors.toCollection(Lists::newArrayList));
//        actualCollapsedReadNames.add(singleFragments.get(0).read().getReadName());
//        Collections.sort(actualCollapsedReadNames);
//        assertEquals(expectedCollapsedReadNames, actualCollapsedReadNames);
//
//        // now check how consensus is formed
//        TestBamWriter writer = new TestBamWriter(config);
//        PartitionReader partitionReader = createPartitionRead(config, writer);
//        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));
//
//        partitionReader.processRead(read1);
//        partitionReader.processRead(read2);
//        partitionReader.processRead(read3);
//        partitionReader.processRead(read4);
//        partitionReader.processRead(read5);
//        partitionReader.processRead(read6);
//
//        partitionReader.postProcessRegion();
//
//        assertEquals(6, writer.nonConsensusWriteCount());
//        assertEquals(1, writer.consensusWriteCount());
//
//        List<SAMRecord> consensusReads = writer.WrittenRecords.stream().filter(x -> x.hasAttribute(CONSENSUS_READ_ATTRIBUTE)).toList();
//        assertEquals(1, consensusReads.size());
//        SAMRecord consensusRead = consensusReads.get(0);
//
//        assertEquals(readName1.split(":")[0], consensusRead.getReadName().split(":")[0]);
//    }

    private String nextReadId(final String umiId)
    {
        return format("%s:%s", mReadIdGen.nextId(), umiId);
    }
}
