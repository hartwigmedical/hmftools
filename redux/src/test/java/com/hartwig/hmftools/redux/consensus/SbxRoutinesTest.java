package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MATCH_SCORE;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MISMATCH_PENALTY;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndels;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.TestUtils.READ_UNMAPPER_DISABLED;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createFragmentCoords;
import static com.hartwig.hmftools.redux.TestUtils.createPartitionRead;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.redux.common.Constants.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.consensus.ConsensusReads.formConsensusReadId;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.getAnnotatedBases;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.processAnnotatedBases;
import static com.hartwig.hmftools.redux.consensus.TemplateReads.selectTemplateRead;
import static com.hartwig.hmftools.redux.umi.UmiConfig.extractUmiIdFromReadId;
import static com.hartwig.hmftools.redux.umi.UmiGroupBuilder.buildUmiGroups;
import static com.hartwig.hmftools.redux.umi.UmiGroupBuilder.collapsePolyGDuplexUmis;
import static com.hartwig.hmftools.redux.umi.UmiGroupBuilder.hasDuplexUmiMatch;
import static com.hartwig.hmftools.redux.umi.UmiUtils.exceedsUmiIdDiff;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.PartitionReader;
import com.hartwig.hmftools.redux.ReadCache;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.TestBamWriter;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.common.ReadInfo;
import com.hartwig.hmftools.redux.umi.PositionFragmentCounts;
import com.hartwig.hmftools.redux.umi.UmiConfig;
import com.hartwig.hmftools.redux.umi.UmiGroupBuilder;
import com.hartwig.hmftools.redux.umi.UmiStatistics;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.jupiter.api.Test;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SbxRoutinesTest
{
    private static final String ZERO_QUAL = String.valueOf(phredToFastq(0));
    private static final String NON_ZERO_QUAL = String.valueOf(phredToFastq(DUPLEX_QUAL));

    @Test
    public void testGetAnnotatedBasesSimple()
    {
        int readLength = 100;
        int alignmentStart = 25;
        String readStr = "A".repeat(readLength);
        String cigar = readLength + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        List<Boolean> duplexIndels = IntStream.range(0, readLength).mapToObj(i -> false).collect(Collectors.toList());

        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);

            assertEquals(i, annotatedBase.ReadIndex);
            assertEquals(i + alignmentStart, annotatedBase.RefPos);
            assertEquals(M, annotatedBase.Op);
            assertEquals('A', (char) annotatedBase.ReadBase);
            assertFalse(annotatedBase.IsDuplexIndel);
        }
    }

    @Test
    public void testGetAnnotatedBasesLeftSoftClip()
    {
        int readLength = 100;
        int alignmentStart = 25;
        int leftSoftClipLength = 10;
        String readStr = "A".repeat(readLength);
        String cigar = leftSoftClipLength + "S" + (readLength - leftSoftClipLength) + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        List<Boolean> duplexIndels = IntStream.range(0, readLength).mapToObj(i -> false).collect(Collectors.toList());

        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);

            CigarOperator expectedOp = i < leftSoftClipLength ? S : M;

            assertEquals(i, annotatedBase.ReadIndex);
            assertEquals(i + alignmentStart - leftSoftClipLength, annotatedBase.RefPos);
            assertEquals(expectedOp, annotatedBase.Op);
            assertEquals('A', (char) annotatedBase.ReadBase);
            assertFalse(annotatedBase.IsDuplexIndel);
        }
    }

    @Test
    public void testGetAnnotatedBasesLeftHardClip()
    {
        int readLength = 100;
        int alignmentStart = 25;
        int leftHardClipLength = 10;
        String readStr = "A".repeat(readLength);
        String cigar = leftHardClipLength + "H" + readLength + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        List<Boolean> duplexIndels = IntStream.range(0, readLength).mapToObj(i -> false).collect(Collectors.toList());

        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);

            assertEquals(i, annotatedBase.ReadIndex);
            assertEquals(i + alignmentStart, annotatedBase.RefPos);
            assertEquals(M, annotatedBase.Op);
            assertEquals('A', (char) annotatedBase.ReadBase);
            assertFalse(annotatedBase.IsDuplexIndel);
        }
    }

    @Test
    public void testGetAnnotatedBasesIndels()
    {
        int readLength = 100;
        int alignmentStart = 25;
        String readStr = "A".repeat(readLength);
        String cigar = "9M1I40M1D50M";
        SAMRecord read = createSamRecordUnpaired(
                "READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);

        List<Boolean> duplexIndels = IntStream.range(0, readLength).mapToObj(i -> false).collect(Collectors.toList());

        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength + 1, annotatedBases.size());

        int expectedReadIndex = 0;
        int expectedRefPos = alignmentStart;
        for(int i = 0; i < annotatedBases.size(); i++)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);

            byte expectedReadBase = i != 50 ? (byte) 'A' : INVALID_BASE_QUAL;
            CigarOperator expectedOp = M;
            if(i == 9)
            {
                expectedOp = I;
            }
            else if(i == 50)
            {
                expectedOp = D;
            }

            assertEquals(expectedReadIndex, annotatedBase.ReadIndex);
            assertEquals(expectedRefPos, annotatedBase.RefPos);
            assertEquals(expectedOp, annotatedBase.Op);
            assertEquals(expectedReadBase, annotatedBase.ReadBase);
            assertFalse(annotatedBase.IsDuplexIndel);

            if(i != 49)
            {
                expectedReadIndex++;
            }

            if(i != 8)
            {
                expectedRefPos++;
            }
        }
    }

    @Test
    public void testProcessAnnotatedBasesForwardRead()
    {
        int alignmentStart = 25;

        String readStr = "A".repeat(5) + "CT".repeat(3) + "A".repeat(5);
        String cigar = "5M4I7M";  // inserts are left aligned.
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);

        String qualStr = NON_ZERO_QUAL.repeat(9) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);
        read.setBaseQualityString(qualStr);

        List<Boolean> duplexIndels = getDuplexIndels("0-9ZZ5-0");

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(1000);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        List<SbxAnnotatedBase> annotedBases = getAnnotatedBases(read, duplexIndels);
        boolean readModified = processAnnotatedBases(refGenome, CHR_1, annotedBases, true);

        List<SbxAnnotatedBase> expectedBases = getAnnotatedBases(read, duplexIndels);

        expectedBases.get(5).setQual((byte) 0);
        expectedBases.get(6).setQual((byte) 0);

        expectedBases.get(5).deleteBase();
        expectedBases.get(6).deleteBase();

        expectedBases.get(7).setQual((byte) 0);
        expectedBases.get(8).setQual((byte) 0);

        assertTrue(readModified);
        assertEquals(expectedBases, annotedBases);
    }

    @Test
    public void testProcessAnnotatedBasesReverseRead()
    {
        int alignmentStart = 25;

        String readStr = "A".repeat(5) + "CT".repeat(3) + "A".repeat(5);
        String cigar = "5M4I7M";  // inserts are left aligned.
        SAMRecord read = createSamRecordUnpaired(
                "READ_001", CHR_1, alignmentStart, readStr, cigar, true, false, null);

        String qualStr = NON_ZERO_QUAL.repeat(5) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(9);
        read.setBaseQualityString(qualStr);

        List<Boolean> duplexIndels = getDuplexIndels("0-9ZZ5-0");
        Collections.reverse(duplexIndels);

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(1000);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        List<SbxAnnotatedBase> annotedBases = getAnnotatedBases(read, duplexIndels);
        boolean readModified = processAnnotatedBases(refGenome, CHR_1, annotedBases, false);

        List<SbxAnnotatedBase> expectedBases = getAnnotatedBases(read, duplexIndels);

        expectedBases.get(9).setQual((byte) 0);
        expectedBases.get(10).setQual((byte) 0);

        expectedBases.get(5).deleteBase();
        expectedBases.get(6).deleteBase();

        expectedBases.get(7).setQual((byte) 0);
        expectedBases.get(8).setQual((byte) 0);

        assertTrue(readModified);
        assertEquals(expectedBases, annotedBases);
    }

    @Test
    public void testStripDuplexIndelsForwardRead()
    {
        int alignmentStart = 25;

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(1000);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        int mapq = 10;
        int nm = 2;
        int alignmentScore = 0;
        String readStr = "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        String cigar = "5M2I7M";  // inserts are left aligned.
        String qualStr = NON_ZERO_QUAL.repeat(7) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);
        String ycTagStr = "0-7ZZ5-0";

        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(mapq);
        read.setAttribute("YC", ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, nm);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignmentScore);

        SbxRoutines.stripDuplexIndels(refGenome, read);

        int alignmentScoreDiff = BWA_GAP_OPEN_PENALTY + 2 * BWA_GAP_EXTEND_PENALTY;

        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(5) + "CT" + "A".repeat(5);
        String expectedCigar = "12M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(5) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);

        SAMRecord expectedRead =
                createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, expectedReadStr, expectedCigar, false, false, null);

        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(mapq);
        expectedRead.setAttribute("YC", ycTagStr);
        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, expectedNm);
        expectedRead.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, expectedAlignmentScore);

        assertEquals(expectedRead, read);
    }

    @Test
    public void testStripDuplexIndelsReverseRead()
    {
        int alignmentStart = 25;

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(1000);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        int mapq = 10;
        int nm = 2;
        int alignmentScore = 0;
        String readStr = "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        String cigar = "5M2I7M";  // inserts are left aligned.
        String qualStr = NON_ZERO_QUAL.repeat(5) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(7);
        String ycTagStr = "0-7ZZ5-0";

        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, true, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(mapq);
        read.setAttribute("YC", ycTagStr);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, nm);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignmentScore);

        SbxRoutines.stripDuplexIndels(refGenome, read);

        int alignmentScoreDiff = BWA_GAP_OPEN_PENALTY + 2 * BWA_GAP_EXTEND_PENALTY;

        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(5) + "CT" + "A".repeat(5);
        String expectedCigar = "12M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(5) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);

        SAMRecord expectedRead =
                createSamRecordUnpaired(
                        "READ_001", CHR_1, alignmentStart, expectedReadStr, expectedCigar, true, false, null);

        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(mapq);
        expectedRead.setAttribute("YC", ycTagStr);
        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, expectedNm);
        expectedRead.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, expectedAlignmentScore);

        assertEquals(expectedRead, read);
    }

    @Test
    public void testFillQualZeroMismatchesWithRefSingleSNV()
    {
        int alignmentStart = 25;
        int mapq = 10;
        int nm = 1;
        int alignmentScore = 0;
        String readStr = "A".repeat(10) + "C" + "A".repeat(9);
        String cigar = "20M";
        String qualStr = NON_ZERO_QUAL.repeat(10) + ZERO_QUAL.repeat(1) + NON_ZERO_QUAL.repeat(9);

        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(mapq);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, nm);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignmentScore);

        String refBases = "A".repeat(alignmentStart - 1 + 10) + "G" + "A".repeat(1000);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        SbxRoutines.fillQualZeroMismatchesWithRef(refGenome, read);

        int alignmentScoreDiff = BWA_MISMATCH_PENALTY + BWA_MATCH_SCORE;

        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(10) + "G" + "A".repeat(9);
        String expectedCigar = "20M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(10) + phredToFastq(1) + NON_ZERO_QUAL.repeat(9);

        SAMRecord expectedRead =
                createSamRecordUnpaired(
                        "READ_001", CHR_1, alignmentStart, expectedReadStr, expectedCigar, false, false, null);
        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(mapq);
        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, expectedNm);
        expectedRead.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, expectedAlignmentScore);

        assertEquals(expectedRead, read);
    }

    private static RefGenomeInterface getRefGenome(final String refBases)
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, refBases);
        mockRefGenome.ChromosomeLengths.put(CHR_1, refBases.length());
        return new CachedRefGenome(mockRefGenome);
    }

    public static class UmiDuplicatesTest
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
            assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(2).DualStrandFrequency);
            assertEquals(1, mPartitionReaderDuplexUMIs.statistics().DuplicateFrequencies.get(3).DualStrandFrequency);
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

            SAMRecord read1 = createSamRecord(
                    "READ_001:" + umiId1, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 200, false, false, null, true, TEST_READ_CIGAR);
            SAMRecord mate1 = createSamRecord(
                    "READ_001:" + umiId1, CHR_2, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200, true, false, null, false, TEST_READ_CIGAR);

            SAMRecord read2 = createSamRecord(
                    "READ_002:" + umiId2, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 200, false, false, null, true, TEST_READ_CIGAR);
            SAMRecord mate2 = createSamRecord(
                    "READ_002:" + umiId2, CHR_2, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200, true, false, null, false, TEST_READ_CIGAR);

            DuplicateGroup readGroup = new DuplicateGroup(Lists.newArrayList(read1, read2), FragmentCoords.fromRead(read1, true));
            DuplicateGroup mateGroup = new DuplicateGroup(Lists.newArrayList(mate2, mate1), FragmentCoords.fromRead(mate1, true));

            List<DuplicateGroup> readUmiGroups = buildUmiGroups(readGroup.fragmentCoordinates(), readGroup.reads(), umiConfig);
            List<DuplicateGroup> mateUmiGroups = buildUmiGroups(mateGroup.fragmentCoordinates(), mateGroup.reads(), umiConfig);

            assertEquals(1, readUmiGroups.size());
            assertEquals(1, mateUmiGroups.size());

            DuplicateGroup readUmiGroup = readUmiGroups.get(0);
            DuplicateGroup mateUmiGroup = mateUmiGroups.get(0);

            SAMRecord readTemplate = TemplateReads.selectTemplateRead(readUmiGroup.reads(), readUmiGroup.fragmentCoordinates());
            String readConsensusReadName = formConsensusReadId(readTemplate, readUmiGroup.umiId());

            SAMRecord mateTemplate = TemplateReads.selectTemplateRead(mateUmiGroup.reads(), mateUmiGroup.fragmentCoordinates());
            String mateConsensusReadName = formConsensusReadId(mateTemplate, mateUmiGroup.umiId());

            assertEquals(readConsensusReadName, mateConsensusReadName);
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

            String umi1 = "AAAAA";
            String umi2 = "AAATA";

            SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                    nextReadId(umi1), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                    nextReadId(umi2), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
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
                    nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
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

            String umidId1Part1 = "TATTAT";
            String umidId2Part1 = "TATGAT";
            String umidIdPart2 = "GCGGCG";
            String umiId = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart2;
            String umiIdReversed = umidIdPart2 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part1;

            SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                    nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
            read1.setFirstOfPairFlag(true);
            read1.setSecondOfPairFlag(false);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(nextReadId(umiIdReversed), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE,
                    TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
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
                    nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
            read1.setFirstOfPairFlag(true);
            read1.setSecondOfPairFlag(false);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(nextReadId(umiIdReversed), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                    1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false, false, null, true, TEST_READ_CIGAR);
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
                    nextReadId(umi), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                    nextReadId(umi), CHR_1, 101, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_001, false, false, null, true, TEST_READ_CIGAR);

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
                    nextReadId(umiId), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
            read1.setFirstOfPairFlag(true);
            read1.setSecondOfPairFlag(false);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(nextReadId(umiIdReversed), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE,
                    TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);
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

            UmiGroupBuilder umiGroupBuilder = new UmiGroupBuilder(ILLUMINA, umiConfig, new UmiStatistics());

            String umi1 = "AAAAA";
            String umi2 = "AAATA";

            SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                    nextReadId(umi1), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                    nextReadId(umi2), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false,
                    false, null, true, TEST_READ_CIGAR);

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


            SAMRecord read1 = SamRecordTestUtils.createSamRecord("READ_001:" + forwardUmiId, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 200, false, false, null, true, TEST_READ_CIGAR);
            read1.setFirstOfPairFlag(true);
            read1.setSecondOfPairFlag(false);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord("READ_002:" + reverseUmiId, CHR_1, 201, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 200, false, false, null, true, TEST_READ_CIGAR);
            read2.setFirstOfPairFlag(false);
            read2.setSecondOfPairFlag(true);

            List<SAMRecord> chr1Reads = Lists.newArrayList(read1, read2);

            SAMRecord mate1 = SamRecordTestUtils.createSamRecord("READ_001:" + forwardUmiId, CHR_2, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200, true, false, null, false, TEST_READ_CIGAR);
            mate1.setFirstOfPairFlag(false);
            mate1.setSecondOfPairFlag(true);

            SAMRecord mate2 = SamRecordTestUtils.createSamRecord("READ_002:" + reverseUmiId, CHR_2, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201, true, false, null, false, TEST_READ_CIGAR);
            mate2.setFirstOfPairFlag(true);
            mate2.setSecondOfPairFlag(false);

            List<SAMRecord> chr2Reads = Lists.newArrayList(mate1, mate2);

            Collections.sort(chr1Reads, Comparator.comparingInt(SAMRecord::getAlignmentStart));
            Collections.sort(chr2Reads, Comparator.comparingInt(SAMRecord::getAlignmentStart));

            // duplicate group forming
            ReadCache readCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, true, ILLUMINA);

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
            UmiGroupBuilder umiGroupBuilder = new UmiGroupBuilder(ILLUMINA, umiConfig, new UmiStatistics());

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
                            upperPos - TEST_READ_BASES.length() + 1, false, false, null, true, TEST_READ_CIGAR);
                    SAMRecord upperRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, upperPos - TEST_READ_BASES.length() + 1,
                            TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, lowerPos, true, false, null, false, TEST_READ_CIGAR);
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
                            upperPos - TEST_READ_BASES.length() + 1, false, false, null, true, TEST_READ_CIGAR);
                    SAMRecord upperRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, upperPos - TEST_READ_BASES.length() + 1,
                            TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, lowerPos, true, false, null, false, TEST_READ_CIGAR);
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
                            upperPos - TEST_READ_BASES.length() + 1, false, false, null, true, TEST_READ_CIGAR);
                    SAMRecord upperRead = SamRecordTestUtils.createSamRecord(readName, CHR_1, upperPos - TEST_READ_BASES.length() + 1,
                            TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, lowerPos, true, false, null, false, TEST_READ_CIGAR);
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

        @Test
        public void testIlluminaPolyGDuplexUmiGroupCollapse()
        {
            MockRefGenome refGenome = new MockRefGenome(true);
            refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
            refGenome.ChromosomeLengths.put(CHR_1, 1_000);

            ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
            TestBamWriter writer = new TestBamWriter(config);
            PartitionReader partitionReader = createPartitionRead(config, writer);

            partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

            String umidIdPart1 = "TCCTATG";
            String umidId1Part2 = "CGGGGGG";
            String umidId2Part2 = "GGGGGGG";
            String umiId1 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
            String umiId2 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;

            final String readName1 = nextReadId(umiId1);
            final String readName2 = nextReadId(umiId2);

            SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                    readName1, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 100, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                    readName2, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100, false, false, null, false, NO_CIGAR);
            read2.setMateUnmappedFlag(true);

            SAMRecord mate2 = SamRecordTestUtils.createSamRecord(
                    readName2, CHR_1, 100, TEST_READ_BASES, NO_CIGAR, CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
            mate2.setReadUnmappedFlag(true);
            flipFirstInPair(mate2);

            partitionReader.processRead(read1.deepCopy());
            partitionReader.processRead(read2.deepCopy());
            partitionReader.processRead(mate2.deepCopy());
            partitionReader.postProcessRegion();

            assertEquals(3, writer.nonConsensusWriteCount());
            assertEquals(0, writer.consensusWriteCount());

            assertEquals(1L, writer.WrittenRecords.stream().filter(x -> x.getReadName().equals(readName1)).count());
            assertEquals(1L, writer.WrittenRecords.stream().filter(x -> x.getReadName().equals(readName2) && x.getFirstOfPairFlag()).count());
            assertEquals(1L, writer.WrittenRecords.stream().filter(x -> x.getReadName().equals(readName2) && x.getSecondOfPairFlag()).count());

            SAMRecord outRead1 = writer.WrittenRecords.stream().filter(x -> x.getReadName().equals(readName1)).findAny().orElse(null);
            SAMRecord outRead2 = writer.WrittenRecords.stream()
                    .filter(x -> x.getReadName().equals(readName2) && x.getFirstOfPairFlag())
                    .findAny()
                    .orElse(null);
            SAMRecord outMate2 = writer.WrittenRecords.stream()
                    .filter(x -> x.getReadName().equals(readName2) && x.getSecondOfPairFlag())
                    .findAny()
                    .orElse(null);

            outRead1.setAttribute(UMI_ATTRIBUTE, null);
            outRead2.setAttribute(UMI_ATTRIBUTE, null);
            outMate2.setAttribute(UMI_ATTRIBUTE, null);

            assertEquals(read1, outRead1);

            assertTrue(outRead2.getDuplicateReadFlag());
            outRead2.setDuplicateReadFlag(false);
            assertEquals(read2, outRead2);

            assertTrue(outMate2.getDuplicateReadFlag());
            outMate2.setDuplicateReadFlag(false);
            assertEquals(mate2, outMate2);
        }

        @Test
        public void testIlluminaPolyGDuplexUmiGroupNoJitterCollapse()
        {
            MockRefGenome refGenome = new MockRefGenome(true);
            refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
            refGenome.ChromosomeLengths.put(CHR_1, 1_000);

            ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
            TestBamWriter writer = new TestBamWriter(config);
            PartitionReader partitionReader = createPartitionRead(config, writer);

            String umidIdPart1 = "TCCTATG";
            String umidId1Part2 = "CGGGGGG";
            String umidId2Part2 = "GGGGGGG";
            String umiId1 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
            String umiId2 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;

            String readName1 = nextReadId(umiId1);
            String readName2 = nextReadId(umiId2);

            SAMRecord read1 = SamRecordTestUtils.createSamRecord(readName1, CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE,
                    TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 100, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(readName2, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                    CHR_1, 100, false, false, null, false, NO_CIGAR);
            read2.setMateUnmappedFlag(true);

            SAMRecord mate2 = SamRecordTestUtils.createSamRecord(readName2, CHR_1, 100, TEST_READ_BASES, NO_CIGAR,
                    CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
            mate2.setReadUnmappedFlag(true);
            flipFirstInPair(mate2);

            partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));
            partitionReader.processRead(read2);
            partitionReader.processRead(mate2);
            partitionReader.processRead(read1);
            partitionReader.postProcessRegion();

            assertEquals(3, writer.nonConsensusWriteCount());
            assertEquals(0, writer.consensusWriteCount());
        }

        @Test
        public void testIlluminaPolyGDuplexUmiGroupNoCollapseWithoutPolyGTail()
        {
            MockRefGenome refGenome = new MockRefGenome(true);
            refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
            refGenome.ChromosomeLengths.put(CHR_1, 1_000);

            ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
            TestBamWriter writer = new TestBamWriter(config);
            PartitionReader partitionReader = createPartitionRead(config, writer);

            partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1_000));

            String umidIdPart1 = "TCCTATG";
            String umidId1Part2 = "CGGGGGA";
            String umidId2Part2 = "GGGGGGA";
            String umiId1 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId1Part2;
            String umiId2 = umidIdPart1 + DEFAULT_DUPLEX_UMI_DELIM + umidId2Part2;

            String readName1 = nextReadId(umiId1);
            SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                    readName1, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 100, false, false, null, true, TEST_READ_CIGAR);

            String readName2 = nextReadId(umiId2);
            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                    readName2, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100, false, false, null, false, NO_CIGAR);
            read2.setMateUnmappedFlag(true);

            SAMRecord mate2 = SamRecordTestUtils.createSamRecord(
                    readName2, CHR_1, 100, TEST_READ_BASES, NO_CIGAR, CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
            mate2.setReadUnmappedFlag(true);
            flipFirstInPair(mate2);

            partitionReader.processRead(read1);
            partitionReader.processRead(read2);
            partitionReader.processRead(mate2);
            partitionReader.postProcessRegion();

            assertEquals(3, writer.nonConsensusWriteCount());
            assertEquals(0, writer.consensusWriteCount());
        }

        @Test
        public void testIlluminaPolyGDuplexUmiGroupNoCollapseWithUmiPrefixMismatch()
        {
            MockRefGenome refGenome = new MockRefGenome(true);
            refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
            refGenome.ChromosomeLengths.put(CHR_1, 1_000);

            ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
            UmiConfig umiConfig = config.UMIs;

            String umidId1Part1 = "TCCTATA";
            String umidId2Part1 = "TCCTATT";
            String umidIdPart2 = "GGGGGGG";
            String umiId1 = umidId1Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart2;
            String umiId2 = umidId2Part1 + DEFAULT_DUPLEX_UMI_DELIM + umidIdPart2;

            String readName1 = nextReadId(umiId1);
            String readName2 = nextReadId(umiId2);

            SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                    readName1, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 100, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                    readName2, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100, false, false, null, false, NO_CIGAR);
            read2.setMateUnmappedFlag(true);

            SAMRecord mate2 = SamRecordTestUtils.createSamRecord(
                    readName2, CHR_1, 100, TEST_READ_BASES, NO_CIGAR, CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
            mate2.setReadUnmappedFlag(true);
            flipFirstInPair(mate2);

            Multiset<SAMRecord> singleFragmentsBefore = HashMultiset.create(List.of(read1.deepCopy(), read2.deepCopy(), mate2.deepCopy()));

            List<DuplicateGroup> umiGroups = Lists.newArrayList();
            List<ReadInfo> singleFragments = Lists.newArrayList(
                    new ReadInfo(read1, FragmentCoords.fromRead(read1, true)),
                    new ReadInfo(read2, FragmentCoords.fromRead(read2, true)),
                    new ReadInfo(mate2, FragmentCoords.fromRead(mate2, true))
            );

            collapsePolyGDuplexUmis(ILLUMINA, umiConfig, umiGroups, singleFragments);
            Multiset<SAMRecord> singleFragmentsAfter = singleFragments.stream()
                    .map(ReadInfo::read)
                    .collect(Collectors.toCollection(HashMultiset::create));

            assertTrue(umiGroups.isEmpty());
            assertEquals(singleFragmentsBefore, singleFragmentsAfter);
        }

        @Test
        public void testIlluminaPolyGDuplexUmiGroupDoNotCollapseTwoFullMappedGroups()
        {
            MockRefGenome refGenome = new MockRefGenome(true);
            refGenome.RefGenomeMap.put(CHR_1, "A".repeat(1_000));
            refGenome.ChromosomeLengths.put(CHR_1, 1_000);

            ReduxConfig config = new ReduxConfig(refGenome, true, true, false, READ_UNMAPPER_DISABLED);
            UmiConfig umiConfig = config.UMIs;

            String umiId1 = "TCCTATG" + DEFAULT_DUPLEX_UMI_DELIM + "GGGGGGG";
            String umiId2 = "TCCTATG" + DEFAULT_DUPLEX_UMI_DELIM + "AAGGGGG";
            String umiId4 = "TCCTATG" + DEFAULT_DUPLEX_UMI_DELIM + "GGGGGGG";

            String readName1 = nextReadId(umiId1);
            String readName2 = nextReadId(umiId2);
            String readName3 = nextReadId(umiId2);
            String readName4 = nextReadId(umiId4);

            SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                    readName1, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 100, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                    readName2, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 100, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read3 = SamRecordTestUtils.createSamRecord(
                    readName3, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 100, false, false, null, true, TEST_READ_CIGAR);

            SAMRecord read4 = SamRecordTestUtils.createSamRecord(
                    readName4, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100, false, false, null, false, NO_CIGAR);
            read4.setMateUnmappedFlag(true);

            SAMRecord mate4 = SamRecordTestUtils.createSamRecord(
                    readName4, CHR_1, 100, TEST_READ_BASES, NO_CIGAR, CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
            mate4.setReadUnmappedFlag(true);
            flipFirstInPair(mate4);

            List<DuplicateGroup> umiGroups = Lists.newArrayList(
                    new DuplicateGroup(umiId1, read1, FragmentCoords.fromRead(read1, true)),
                    new DuplicateGroup(umiId2, Lists.newArrayList(read2, read3), FragmentCoords.fromRead(read2, true))
            );
            List<ReadInfo> singleFragments = Lists.newArrayList(
                    new ReadInfo(read4, FragmentCoords.fromRead(read4, true)),
                    new ReadInfo(mate4, FragmentCoords.fromRead(mate4, true))
            );

            collapsePolyGDuplexUmis(ILLUMINA, umiConfig, umiGroups, singleFragments);

            assertEquals(1, umiGroups.size());

            DuplicateGroup umiGroup = umiGroups.get(0);

            assertEquals(HashMultiset.create(List.of(readName2, readName3)),
                    umiGroup.reads().stream().map(SAMRecord::getReadName).collect(Collectors.toCollection(HashMultiset::create)));
            assertEquals(HashMultiset.create(List.of(readName4, readName4)),
                    umiGroup.polyGUmiReads().stream().map(SAMRecord::getReadName).collect(Collectors.toCollection(HashMultiset::create)));

            assertEquals(HashMultiset.create(List.of(readName1)),
                    singleFragments.stream()
                            .map(ReadInfo::read)
                            .map(SAMRecord::getReadName)
                            .collect(Collectors.toCollection(HashMultiset::create)));
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

            SAMRecord read1 = createSamRecord(readName1, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 1_000, false, false, null, false, TEST_READ_CIGAR);
            flipFirstInPair(read1);

            SAMRecord mate1 = createSamRecord(readName1, CHR_2, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, false, TEST_READ_CIGAR);

            SAMRecord read2 = createSamRecord(readName2, CHR_2, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 1_000, false, false, null, false, NO_CIGAR);
            read2.setMateUnmappedFlag(true);

            SAMRecord mate2 = createSamRecord(readName2, CHR_2, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_2, 1_000, false, false, null, false, TEST_READ_CIGAR);
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
                    readStream = (final String readName, final Boolean isFirstOfPair) -> writtenRecords.stream().filter(x -> (x.getReadName().equals(readName) && x.getFirstOfPairFlag() == isFirstOfPair));

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
                    readName1, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false, false, null, false, TEST_READ_CIGAR);
            read1.setReadUnmappedFlag(true);

            SAMRecord mate1 = createSamRecord(
                    readName1, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, false, NO_CIGAR);
            mate1.setMateUnmappedFlag(true);
            flipFirstInPair(mate1);

            SAMRecord read2 = createSamRecord(
                    readName2, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false, false, null, false, TEST_READ_CIGAR);
            read2.setReadUnmappedFlag(true);

            SAMRecord mate2 = createSamRecord(
                    readName2, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, false, NO_CIGAR);
            mate2.setMateUnmappedFlag(true);
            flipFirstInPair(mate2);

            SAMRecord read3 = createSamRecord(
                    readName3, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false, false, null, false, TEST_READ_CIGAR);
            read3.setReadUnmappedFlag(true);

            SAMRecord mate3 = createSamRecord(
                    readName3, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, false, NO_CIGAR);
            mate3.setMateUnmappedFlag(true);
            flipFirstInPair(mate3);

            SAMRecord read4 = createSamRecord(
                    readName4, CHR_1, 1_000, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false, false, null, false, TEST_READ_CIGAR);
            read4.setReadUnmappedFlag(true);

            SAMRecord mate4 = createSamRecord(
                    readName1, CHR_1, 1_000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_000, false, false, null, false, NO_CIGAR);
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
                    readName1, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false, false, null, false, NO_CIGAR);
            read1.setMateUnmappedFlag(true);

            SAMRecord mate1 = createSamRecord(
                    readName1, CHR_1, 500, TEST_READ_BASES, NO_CIGAR, CHR_1, 500, false, false, null, false, TEST_READ_CIGAR);
            mate1.setReadUnmappedFlag(true);
            flipFirstInPair(mate1);

            SAMRecord read2 = createSamRecord(
                    readName2, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false, false, null, false, NO_CIGAR);
            read2.setMateUnmappedFlag(true);

            SAMRecord mate2 = createSamRecord(
                    readName2, CHR_1, 500, TEST_READ_BASES, NO_CIGAR, CHR_1, 500, false, false, null, false, TEST_READ_CIGAR);
            mate2.setReadUnmappedFlag(true);
            flipFirstInPair(mate2);

            SAMRecord read3 = createSamRecord(
                    readName3, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_500, false, false, null, false, TEST_READ_CIGAR);
            SAMRecord mate3 = createSamRecord(
                    readName3, CHR_1, 1_500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false, false, null, false, TEST_READ_CIGAR);
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
                    readName4, CHR_1, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1_500, false, false, null, false, TEST_READ_CIGAR);
            SAMRecord mate4 = createSamRecord(
                    readName4, CHR_1, 1_500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 500, false, false, null, false, TEST_READ_CIGAR);
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

        private String nextReadId(final String umiId)
        {
            return format("%s:%s", mReadIdGen.nextId(), umiId);
        }
    }

    public static class UmiGroupsTest
    {
        public static final String FIXED_READ_ID = "123:ABC:1:4455:";

        private static final UmiConfig UMI_CONFIG = new UmiConfig(true, false, "", false);

        @Test
        public void testUmiUtils()
        {
            String umiId = "TATCGC";
            String readId = FIXED_READ_ID + umiId;
            assertEquals(umiId, extractUmiIdFromReadId(readId));

            short umiLength = (short)umiId.length();
            assertEquals(umiId, UMI_CONFIG.extractUmiId(readId, umiLength));

            String readId1 = FIXED_READ_ID + "TATCGC";
            String readId2 = FIXED_READ_ID + "TATCGG";
            String readId3 = FIXED_READ_ID + "TATGGG";
            String readId4 = FIXED_READ_ID + "AAACGG";
            String readId5 = FIXED_READ_ID + "AAAGGG";
            String readId6 = FIXED_READ_ID + "AAAGGGT";

            assertFalse(exceedsUmiIdDiff(readId1, readId2));
            assertFalse(exceedsUmiIdDiff(readId2, readId3));
            assertTrue(exceedsUmiIdDiff(readId1, readId3));
            assertTrue(exceedsUmiIdDiff(readId2, readId4));
            assertFalse(exceedsUmiIdDiff(readId4, readId5));
            assertTrue(exceedsUmiIdDiff(readId5, readId6));
        }

        @Test
        public void testDuplexUmis()
        {
            UmiConfig umiConfig = new UmiConfig(true, true, String.valueOf(DEFAULT_DUPLEX_UMI_DELIM), false);

            String umiId1 = "TATCGC_AAGTCG";
            assertFalse(hasDuplexUmiMatch(umiId1, umiId1, umiConfig.DuplexDelim, umiConfig.PermittedBaseDiff));

            String umiId2 = "AAGTCC_TATCGG";

            assertTrue(hasDuplexUmiMatch(umiId1, umiId2, umiConfig.DuplexDelim, umiConfig.PermittedBaseDiff));

            // too manu different base
            umiId2 = "AAGCCC_TATCGC";

            assertFalse(hasDuplexUmiMatch(umiId1, umiId2, umiConfig.DuplexDelim, umiConfig.PermittedBaseDiff));
        }

        private static SAMRecord createSimpleRead(final String readId)
        {
            return createSamRecord(readId, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                    false, false, null, true, TEST_READ_CIGAR);
        }

        @Test
        public void testUmiGroupAssignment()
        {
            SAMRecord frag1 = createSimpleRead(FIXED_READ_ID + "TATCGC");
            SAMRecord frag2 = createSimpleRead(FIXED_READ_ID + "TATCGC");
            SAMRecord frag3 = createSimpleRead(FIXED_READ_ID + "TATCCC");
            SAMRecord frag4 = createSimpleRead(FIXED_READ_ID + "AATGGG");
            SAMRecord frag5 = createSimpleRead(FIXED_READ_ID + "TATCGC");
            SAMRecord frag6 = createSimpleRead(FIXED_READ_ID + "AATGGG");

            List<SAMRecord> fragments = Lists.newArrayList(frag1, frag2, frag3, frag4, frag5, frag6);
            FragmentCoords fragmentCoords = createFragmentCoords(frag1);
            List<DuplicateGroup> groups = buildUmiGroups(fragmentCoords, fragments, UMI_CONFIG);
            assertEquals(2, groups.size());

            DuplicateGroup group = groups.stream().filter(x -> x.reads().contains(frag1)).findFirst().orElse(null);
            assertTrue(group.reads().contains(frag2));
            assertTrue(group.reads().contains(frag3));
            assertTrue(group.reads().contains(frag5));

            group = groups.stream().filter(x -> x.reads().contains(frag4)).findFirst().orElse(null);
            assertTrue(group.reads().contains(frag6));

            // check that larger groups aren't merged
            SAMRecord frag11 = createSimpleRead(FIXED_READ_ID + "TTTCGT");
            SAMRecord frag12 = createSimpleRead(FIXED_READ_ID + "TTCCGT");
            SAMRecord frag13 = createSimpleRead(FIXED_READ_ID + "TTACGT");
            SAMRecord frag14 = createSimpleRead(FIXED_READ_ID + "TTACAT");
            SAMRecord frag15 = createSimpleRead(FIXED_READ_ID + "TAAAAT");
            SAMRecord frag16 = createSimpleRead(FIXED_READ_ID + "TTACAG");

            // duplicate some to build bigger UMI groups
            fragments = Lists.newArrayList(
                    frag11, frag11,
                    frag12, frag12, frag12,
                    frag13, frag13, frag13,frag13, frag13, frag13, frag13, frag13,
                    frag14, frag14, frag14, frag14,
                    frag15, frag15, frag15, frag15, frag15,
                    frag16);

            groups = buildUmiGroups(fragmentCoords, fragments, UMI_CONFIG);
            assertEquals(2, groups.size());
            assertEquals(fragments.size(), groups.stream().mapToInt(x -> x.readCount()).sum());

            group = groups.stream().filter(x -> x.reads().contains(frag11)).findFirst().orElse(null);
            assertEquals(18, group.readCount());
            assertTrue(group.reads().contains(frag12));
            assertTrue(group.reads().contains(frag13));
            assertTrue(group.reads().contains(frag14));
            assertTrue(group.reads().contains(frag16));

            group = groups.stream().filter(x -> x.reads().contains(frag15)).findFirst().orElse(null);
            assertEquals(5, group.readCount());
        }

        @Test
        public void testUmiGroupAssignment2()
        {
            // test a final merge of groups with 2 bases difference
            SAMRecord frag1 = createSimpleRead(FIXED_READ_ID + "TTAAGG");
            SAMRecord frag2 = createSimpleRead(FIXED_READ_ID + "TTAAGC");
            SAMRecord frag3 = createSimpleRead(FIXED_READ_ID + "TTAATT");
            SAMRecord frag4 = createSimpleRead(FIXED_READ_ID + "GGGATT");
            SAMRecord frag5 = createSimpleRead(FIXED_READ_ID + "GGGAGG");
            SAMRecord frag6 = createSimpleRead(FIXED_READ_ID + "CCGAGC");

            List<SAMRecord> fragments = Lists.newArrayList(frag1, frag2, frag3, frag4, frag5, frag6);
            FragmentCoords fragmentCoords = createFragmentCoords(frag1);
            List<DuplicateGroup> groups = buildUmiGroups(fragmentCoords, fragments, UMI_CONFIG);
            assertEquals(3, groups.size());
        }

        @Test
        public void testUmiGroupAssignment3()
        {
            // test a final merge of groups with 4 bases difference if the groups are large enough
            String umi1 = "TTAAGG";
            String umi2 = "TTAGGA"; // 2-base diff from #1
            String umi3 = "CCGGTT"; // unrelated
            String umi4 = "TTCCCC"; // 4-base diff from #1 and #2

            List<SAMRecord> fragments = Lists.newArrayList();

            fragments.add(createSimpleRead(FIXED_READ_ID + umi1));
            fragments.add(createSimpleRead(FIXED_READ_ID + umi2));
            fragments.add(createSimpleRead(FIXED_READ_ID + umi3));

            for(int i = 0; i < 51; ++i)
            {
                fragments.add(createSimpleRead(FIXED_READ_ID + umi4));
            }

            FragmentCoords fragmentCoords = createFragmentCoords(fragments.get(0));
            List<DuplicateGroup> groups = buildUmiGroups(fragmentCoords, fragments, UMI_CONFIG);
            assertEquals(2, groups.size());
        }

        @Test
        public void testDefinedUmis()
        {
            UmiConfig umiConfig = new UmiConfig(true, false, "", false);
            String definedUmi1 = "AAAGGG";
            String definedUmi2 = "TTTAAA";
            String definedUmi3 = "CCCAAA";
            umiConfig.addDefinedUmis(Sets.newHashSet(definedUmi1, definedUmi2, definedUmi3));

            SAMRecord frag1 = createSimpleRead(FIXED_READ_ID + definedUmi1);
            SAMRecord frag2 = createSimpleRead(FIXED_READ_ID + "AAAGGC");
            SAMRecord frag3 = createSimpleRead(FIXED_READ_ID + definedUmi1);
            SAMRecord frag4 = createSimpleRead(FIXED_READ_ID + definedUmi2);
            SAMRecord frag5 = createSimpleRead(FIXED_READ_ID + "TTTAAC");
            SAMRecord frag6 = createSimpleRead(FIXED_READ_ID + definedUmi3);

            List<SAMRecord> fragments = Lists.newArrayList(frag1, frag2, frag3, frag4, frag5, frag6);
            FragmentCoords fragmentCoords = createFragmentCoords(fragments.get(0));
            List<DuplicateGroup> groups = buildUmiGroups(fragmentCoords, fragments, umiConfig);
            assertEquals(3, groups.size());
        }

        /*
        @Test
        public void testPerfUmiIdExtraction()
        {
            String umiId = "TATCGC";
            String readId = FIXED_READ_ID + umiId;
            short umiLength = (short) umiId.length();

            int testIterations = 1000000;
            long startTimeMs = System.currentTimeMillis();

            for(int i = 0; i < testIterations; ++i)
            {
                extractUmiId(readId);
            }

            long timeTakenMs = System.currentTimeMillis() - startTimeMs;

            BM_LOGGER.info("slow time: {}", timeTakenMs);

            startTimeMs = System.currentTimeMillis();

            for(int i = 0; i < testIterations; ++i)
            {
                extractUmiId(readId, umiLength);
            }

            long timeTakenMs2 = System.currentTimeMillis() - startTimeMs;

            BM_LOGGER.info("fast time: {}", timeTakenMs2);
        }
        */
    }
}
