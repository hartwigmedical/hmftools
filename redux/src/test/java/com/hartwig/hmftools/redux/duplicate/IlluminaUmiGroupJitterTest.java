package com.hartwig.hmftools.redux.duplicate;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.ReduxConstants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.redux.ReduxConstants.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.TestUtils.READ_UNMAPPER_DISABLED;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createPartitionRead;
import static com.hartwig.hmftools.redux.consensus.TemplateReads.selectTemplateRead;
import static com.hartwig.hmftools.redux.duplicate.UmiDuplicatesTest.nextReadId;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.PartitionReader;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.TestBamWriter;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class IlluminaUmiGroupJitterTest
{
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

        String umi1 = "AAAAA";
        String umi2 = "AAATA";

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi1), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi2), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                1_000, false,
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

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false,
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

        String umidId1Part1 = "TATTAT";
        String umidId2Part1 = "TATGAT";
        String umidIdPart2 = "GCGGCG";
        String umiId = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart2;
        String umiIdReversed = umidIdPart2 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part1;

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(nextReadId(umiIdReversed), CHR_1,
                100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE,
                TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null,
                true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(false);
        read1.setSecondOfPairFlag(true);

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

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(nextReadId(umiIdReversed), CHR_1, 100, TEST_READ_BASES,
                TEST_READ_CIGAR, CHR_1, 1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false, false,
                null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(false);
        read1.setSecondOfPairFlag(true);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
        assertEquals(1, writer.consensusWriteCount());
    }

    @Test
    public void testIlluminaJitterUmiGroupNoCollapseDueToNoFixedFragmentEnd()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, false, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

        String umi = "AAAAA";

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);
        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi), CHR_1, 101, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_001, false,
                false, null, true, TEST_READ_CIGAR);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
        assertEquals(0, writer.consensusWriteCount());
    }

    @Test
    public void testIlluminaJitterDuplexUmiGroupNoCollapseDueToUMIMismatch()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

        String umidIdPart1 = "TATTAT";
        String umidId1Part2 = "GCGGCG";
        String umidId2Part2 = "GCATCG";
        String umiId = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
        String umiIdReversed = umidId2Part2 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart1;

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(nextReadId(umiIdReversed), CHR_1,
                100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000,
                false, false, null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(false);
        read1.setSecondOfPairFlag(true);

        partitionReader.processRead(read1);
        partitionReader.processRead(read2);
        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
        assertEquals(0, writer.consensusWriteCount());
    }

    @Test
    public void testIlluminaJitterUmiGroupCollapsedReadsNotUsedForConsensus()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, false, false, READ_UNMAPPER_DISABLED);
        UmiConfig umiConfig = config.UMIs;

        UmiGroupBuilder umiGroupBuilder = new UmiGroupBuilder(umiConfig, new UmiStatistics());

        String umi1 = "AAAAA";
        String umi2 = "AAATA";

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi1), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                false, null, true, TEST_READ_CIGAR);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                nextReadId(umi2), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                1_000, false, false, null, true, TEST_READ_CIGAR);

        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();
        List<ReadInfo> singleFragments = Lists.newArrayList(
                new ReadInfo(read1, FragmentCoords.fromRead(read1, true)),
                new ReadInfo(read2, FragmentCoords.fromRead(read2, true))
        );

        List<DuplicateGroup> umiGroups = umiGroupBuilder.processUmiGroups(duplicateGroups, singleFragments, true);

        assertTrue(singleFragments.isEmpty());
        assertEquals(1, umiGroups.size());

        DuplicateGroup umiGroup = umiGroups.get(0);

        assertEquals(1, umiGroup.reads().size());
        assertEquals(1, umiGroup.nonConsensusReads().size());
    }

    @Test
    public void testConsistentConsensusReadPairAfterIlluminaJitterUmiGroupCollapse()
    {
        String umiIdPart1 = "A".repeat(6);
        String umiIdPart2 = "C".repeat(6);
        String forwardUmiId = umiIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umiIdPart2;
        String reverseUmiId = umiIdPart2 + DEFAULT_DUPLEX_UMI_DELIM + umiIdPart1;


        SAMRecord read1 = SamRecordTestUtils.createSamRecord("READ_001:" + forwardUmiId, CHR_1, 200, TEST_READ_BASES,
                TEST_READ_CIGAR, CHR_2, 200, false, false, null, true, TEST_READ_CIGAR);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord("READ_002:" + reverseUmiId, CHR_1, 201, TEST_READ_BASES,
                TEST_READ_CIGAR, CHR_2, 200, false, false, null, true, TEST_READ_CIGAR);
        read2.setFirstOfPairFlag(false);
        read2.setSecondOfPairFlag(true);

        List<SAMRecord> chr1Reads = Lists.newArrayList(read1, read2);

        SAMRecord mate1 = SamRecordTestUtils.createSamRecord("READ_001:" + forwardUmiId, CHR_2, 200, TEST_READ_BASES,
                TEST_READ_CIGAR, CHR_1, 200, true, false, null, false, TEST_READ_CIGAR);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        SAMRecord mate2 = SamRecordTestUtils.createSamRecord("READ_002:" + reverseUmiId, CHR_2, 200, TEST_READ_BASES,
                TEST_READ_CIGAR, CHR_1, 201, true, false, null, false, TEST_READ_CIGAR);
        mate2.setFirstOfPairFlag(true);
        mate2.setSecondOfPairFlag(false);

        List<SAMRecord> chr2Reads = Lists.newArrayList(mate1, mate2);

        Collections.sort(chr1Reads, Comparator.comparingInt(SAMRecord::getAlignmentStart));
        Collections.sort(chr2Reads, Comparator.comparingInt(SAMRecord::getAlignmentStart));

        // duplicate group forming
        ReadCache readCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, true);

        chr1Reads.forEach(readCache::processRead);
        FragmentCoordReads chr1FragmentCoordsReads = readCache.evictAll();
        readCache.clear();

        chr2Reads.forEach(readCache::processRead);
        FragmentCoordReads chr2FragmentCoordsReads = readCache.evictAll();
        readCache.clear();

        assertEquals(0, chr1FragmentCoordsReads.DuplicateGroups.size());
        assertEquals(2, chr1FragmentCoordsReads.SingleReads.size());
        assertEquals(2, chr1FragmentCoordsReads.totalReadCount());

        assertEquals(0, chr2FragmentCoordsReads.DuplicateGroups.size());
        assertEquals(2, chr2FragmentCoordsReads.SingleReads.size());
        assertEquals(2, chr2FragmentCoordsReads.totalReadCount());

        // umi group forming
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_1, 1_000);
        refGenome.RefGenomeMap.put(CHR_2, "A".repeat(1_000));
        refGenome.ChromosomeLengths.put(CHR_2, 1_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
        UmiConfig umiConfig = config.UMIs;
        UmiGroupBuilder umiGroupBuilder = new UmiGroupBuilder(umiConfig, new UmiStatistics());

        List<ReadInfo> chr1SingleReads = chr1FragmentCoordsReads.SingleReads;
        List<DuplicateGroup> chr1UmiGroups = umiGroupBuilder.processUmiGroups(
                chr1FragmentCoordsReads.DuplicateGroups, chr1SingleReads, true);

        List<ReadInfo> chr2SingleReads = chr2FragmentCoordsReads.SingleReads;
        List<DuplicateGroup> chr2UmiGroups = umiGroupBuilder.processUmiGroups(
                chr2FragmentCoordsReads.DuplicateGroups, chr2SingleReads, true);

        assertEquals(1, chr1UmiGroups.size());
        assertEquals(0, chr1SingleReads.size());

        assertEquals(1, chr2UmiGroups.size());
        assertEquals(0, chr2SingleReads.size());

        List<String> chr1SingleReadNames = chr1SingleReads.stream().map(x -> x.read().getReadName()).sorted().collect(Collectors.toList());
        List<String> chr2SingleReadNames = chr2SingleReads.stream().map(x -> x.read().getReadName()).sorted().collect(Collectors.toList());
        assertEquals(chr1SingleReadNames, chr2SingleReadNames);

        List<String> chr1TemplateReadNames = chr1UmiGroups.stream()
                .map(x -> selectTemplateRead(x.reads(), x.fragmentCoordinates()).getReadName())
                .sorted()
                .collect(Collectors.toList());

        List<String> chr2TemplateReadNames = chr2UmiGroups.stream()
                .map(x -> selectTemplateRead(x.reads(), x.fragmentCoordinates()).getReadName())
                .sorted()
                .collect(Collectors.toList());

        assertEquals(chr1TemplateReadNames, chr2TemplateReadNames);

        final Function<List<SAMRecord>, Set<String>> readsToReadNameSet =
                xs -> xs.stream().map(SAMRecord::getReadName).collect(Collectors.toCollection(Sets::newHashSet));

        Set<Pair<Set<String>, Set<String>>> chr1UmiGroupsSet = chr1UmiGroups.stream()
                .map(x -> Pair.of(readsToReadNameSet.apply(x.reads()), readsToReadNameSet.apply(x.nonConsensusReads())))
                .collect(Collectors.toCollection(Sets::newHashSet));

        Set<Pair<Set<String>, Set<String>>> chr2UmiGroupsSet = chr2UmiGroups.stream()
                .map(x -> Pair.of(readsToReadNameSet.apply(x.reads()), readsToReadNameSet.apply(x.nonConsensusReads())))
                .collect(Collectors.toCollection(Sets::newHashSet));

        assertEquals(chr1UmiGroupsSet, chr2UmiGroupsSet);
    }

    @Test
    public void testIlluminaJitterUmiGroupCollapseReadCacheInconsistentPops()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(2_000));
        refGenome.ChromosomeLengths.put(CHR_1, 2_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, false, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 2_000));

        String umiId1 = "CATGATC";
        String umiId2 = "TCCTATC";

        record Fragment(String readName, int lowerPos, int upperPos)
        {
            private Collection<SAMRecord> samRecords()
            {
                SAMRecord lowerRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, lowerPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                        upperPos - TEST_READ_BASES.length() + 1, false, false, null,
                        true, TEST_READ_CIGAR);

                SAMRecord upperRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, upperPos - TEST_READ_BASES.length() + 1,
                        TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, lowerPos, true, false, null,
                        false, TEST_READ_CIGAR);

                flipFirstInPair(upperRead);

                return List.of(lowerRead, upperRead);
            }

            public static Collection<SAMRecord> createSamRecords(final String readName, int lowerPos, int upperPos)
            {
                Fragment fragment = new Fragment(readName, lowerPos, upperPos);
                return fragment.samRecords();
            }
        }

        List<SAMRecord> reads = Lists.newArrayList();
        reads.addAll(Fragment.createSamRecords("READ_001:" + umiId1, 100, 700));
        reads.addAll(Fragment.createSamRecords("READ_002:" + umiId1, 101, 700));
        reads.addAll(Fragment.createSamRecords("READ_003:" + umiId2, 100 + ReadCache.DEFAULT_MAX_SOFT_CLIP, 1_000));

        reads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        reads.forEach(read -> partitionReader.processRead(read));

        partitionReader.postProcessRegion();

        assertEquals(6, writer.nonConsensusWriteCount());
        assertEquals(2, writer.consensusWriteCount());

        Multiset<String> readNames = writer.WrittenRecords.stream()
                .map(SAMRecord::getReadName)
                .collect(Collectors.toCollection(HashMultiset::create));

        assertEquals(8, readNames.size());
        assertEquals(0, (int) readNames.entrySet().stream().mapToInt(Multiset.Entry::getCount).filter(x -> x != 2).count());
    }

    @Test
    public void testIlluminaJitterUmiGroupCollapseReadCacheInconsistentPopsNotFixedByIncreasingMaxSoftClip()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(2_000));
        refGenome.ChromosomeLengths.put(CHR_1, 2_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, false, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 2_000));

        String umiId1 = "CATGATC";
        String umiId2 = "TCCTATC";

        record Fragment(String readName, int lowerPos, int upperPos)
        {
            private Collection<SAMRecord> samRecords()
            {
                SAMRecord lowerRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, lowerPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                        upperPos - TEST_READ_BASES.length() + 1, false, false, null,
                        true, TEST_READ_CIGAR);

                SAMRecord upperRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, upperPos - TEST_READ_BASES.length() + 1,
                        TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, lowerPos, true, false, null,
                        false, TEST_READ_CIGAR);
                flipFirstInPair(upperRead);

                return List.of(lowerRead, upperRead);
            }

            public static Collection<SAMRecord> createSamRecords(final String readName, int lowerPos, int upperPos)
            {
                Fragment fragment = new Fragment(readName, lowerPos, upperPos);
                return fragment.samRecords();
            }
        }

        List<SAMRecord> reads = Lists.newArrayList();
        reads.addAll(Fragment.createSamRecords("READ_001:" + umiId1, 100, 700));
        reads.addAll(Fragment.createSamRecords("READ_002:" + umiId1, 101, 700));
        reads.addAll(Fragment.createSamRecords(
                "READ_003:" + umiId2, 100 + ReadCache.DEFAULT_MAX_SOFT_CLIP + SINGLE_END_JITTER_COLLAPSE_DISTANCE, 1_000));

        reads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        reads.forEach(read -> partitionReader.processRead(read));

        partitionReader.postProcessRegion();

        assertEquals(6, writer.nonConsensusWriteCount());
        assertEquals(2, writer.consensusWriteCount());

        Multiset<String> readNames = writer.WrittenRecords.stream()
                .map(SAMRecord::getReadName)
                .collect(Collectors.toCollection(HashMultiset::create));

        assertEquals(8, readNames.size());
        assertEquals(0, (int) readNames.entrySet().stream().mapToInt(Multiset.Entry::getCount).filter(x -> x != 2).count());
    }

    @Test
    public void testIlluminaJitterUmiGroupCollapseReadCacheInconsistentPopsNotFixedByIncreasingMaxSoftClipAsAnInitialFilter()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(2_000));
        refGenome.ChromosomeLengths.put(CHR_1, 2_000);

        ReduxConfig config = new ReduxConfig(refGenome, true, false, false, READ_UNMAPPER_DISABLED);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 2_000));

        String umiId1 = "CATGATC";
        String umiId2 = "TCCTATC";

        record Fragment(String readName, int lowerPos, int upperPos)
        {
            private Collection<SAMRecord> samRecords()
            {
                SAMRecord lowerRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, lowerPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                        upperPos - TEST_READ_BASES.length() + 1, false, false, null,
                        true, TEST_READ_CIGAR);
                SAMRecord upperRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, upperPos - TEST_READ_BASES.length() + 1,
                        TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, lowerPos, true, false, null,
                        false, TEST_READ_CIGAR);
                flipFirstInPair(upperRead);

                return List.of(lowerRead, upperRead);
            }

            public static Collection<SAMRecord> createSamRecords(final String readName, int lowerPos, int upperPos)
            {
                Fragment fragment = new Fragment(readName, lowerPos, upperPos);
                return fragment.samRecords();
            }
        }

        List<SAMRecord> reads = Lists.newArrayList();
        reads.addAll(Fragment.createSamRecords("READ_001:" + umiId1, 100, 700));
        reads.addAll(Fragment.createSamRecords("READ_002:" + umiId1, 110, 700));
        reads.addAll(Fragment.createSamRecords("READ_003:" + umiId1, 111, 700));
        reads.addAll(Fragment.createSamRecords("READ_004:" + umiId1, 111, 700));
        reads.addAll(Fragment.createSamRecords(
                "READ_005:" + umiId2, 100 + ReadCache.DEFAULT_MAX_SOFT_CLIP + SINGLE_END_JITTER_COLLAPSE_DISTANCE, 1_000));

        reads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        reads.forEach(read -> partitionReader.processRead(read));

        partitionReader.postProcessRegion();

        assertEquals(10, writer.nonConsensusWriteCount());
        assertEquals(2, writer.consensusWriteCount());

        Multiset<String> readNames = writer.WrittenRecords.stream()
                .map(SAMRecord::getReadName)
                .collect(Collectors.toCollection(HashMultiset::create));

        assertEquals(12, readNames.size());
        assertEquals(0, (int) readNames.entrySet().stream().mapToInt(Multiset.Entry::getCount).filter(x -> x != 2).count());
    }

}
