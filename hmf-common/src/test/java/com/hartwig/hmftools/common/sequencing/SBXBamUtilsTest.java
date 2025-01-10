package com.hartwig.hmftools.common.sequencing;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.ALIGNMENT_SCORED_DIFF_TO_MAPQ_DIFF;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.BWA_MATCH_SCORE;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.BWA_MISMATCH_PENALTY;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.MAX_MAPQ;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.fillQualZeroMismatchesWithRef;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.getAnnotatedBases;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.getDuplexIndels;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.processAnnotatedBases;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.stripDuplexIndels;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.sequencing.SBXBamUtils.AnnotatedBase;

import org.junit.Test;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SBXBamUtilsTest
{
    private static final String ZERO_QUAL = String.valueOf(phredToFastq(0));
    private static final String NON_ZERO_QUAL = String.valueOf(phredToFastq(DUPLEX_QUAL));

    @Test
    public void testGetDuplexIndelsNoDuplexIndels()
    {
        String ycTagStr = "0-100-0";
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(100, duplexIndels.size());
        assertTrue(duplexIndels.stream().noneMatch(x -> x));
    }

    @Test
    public void testGetDuplexIndelsNoDuplexIndelsWithSimplexHead()
    {
        String ycTagStr = "10-100-0";
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(110, duplexIndels.size());
        assertTrue(duplexIndels.stream().noneMatch(x -> x));
    }

    @Test
    public void testGetDuplexIndelsSingleDuplexIndel()
    {
        String ycTagStr = "0-100I100-0";
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(201, duplexIndels.size());
        assertTrue(duplexIndels.subList(0, 100).stream().noneMatch(x -> x));
        assertTrue(duplexIndels.get(100));
        assertTrue(duplexIndels.subList(101, 201).stream().noneMatch(x -> x));
    }

    @Test
    public void testGetDuplexIndelsSingleNonDuplexIndel()
    {
        String ycTagStr = "0-100A100-0";
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(201, duplexIndels.size());
        assertTrue(duplexIndels.stream().noneMatch(x -> x));
    }

    @Test
    public void testGetAnnotatedBasesSimple()
    {
        int readLength = 100;
        int alignmentStart = 25;
        String readStr = "A".repeat(readLength);
        String cigar = readLength + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        List<Boolean> duplexIndels = IntStream.range(0, readLength).mapToObj(i -> false).collect(Collectors.toList());

        List<AnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            AnnotatedBase annotatedBase = annotatedBases.get(i);

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

        List<AnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            AnnotatedBase annotatedBase = annotatedBases.get(i);

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

        List<AnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength, annotatedBases.size());

        for(int i = 0; i < annotatedBases.size(); i++)
        {
            AnnotatedBase annotatedBase = annotatedBases.get(i);

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
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        List<Boolean> duplexIndels = IntStream.range(0, readLength).mapToObj(i -> false).collect(Collectors.toList());

        List<AnnotatedBase> annotatedBases = getAnnotatedBases(read, duplexIndels);

        assertEquals(readLength + 1, annotatedBases.size());

        int expectedReadIndex = 0;
        int expectedRefPos = alignmentStart;
        for(int i = 0; i < annotatedBases.size(); i++)
        {
            AnnotatedBase annotatedBase = annotatedBases.get(i);

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

        String readStr = "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        String cigar = "5M2I7M";  // inserts are left aligned.
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);

        String qualStr = NON_ZERO_QUAL.repeat(7) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);
        read.setBaseQualityString(qualStr);

        List<Boolean> duplexIndels = getDuplexIndels("0-7ZZ5-0");

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(1000);

        List<AnnotatedBase> annotedBases = getAnnotatedBases(read, duplexIndels);
        boolean readModified = processAnnotatedBases(annotedBases, true, refBases.getBytes(), 1);

        List<AnnotatedBase> expectedBases = getAnnotatedBases(read, duplexIndels);

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

        String readStr = "A".repeat(5) + "CT".repeat(2) + "A".repeat(5);
        String cigar = "5M2I7M";  // inserts are left aligned.
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, true, false, null);

        String qualStr = NON_ZERO_QUAL.repeat(5) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(7);
        read.setBaseQualityString(qualStr);

        List<Boolean> duplexIndels = getDuplexIndels("0-7ZZ5-0");
        Collections.reverse(duplexIndels);

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(1000);

        List<AnnotatedBase> annotedBases = getAnnotatedBases(read, duplexIndels);
        boolean readModified = processAnnotatedBases(annotedBases, false, refBases.getBytes(), 1);

        List<AnnotatedBase> expectedBases = getAnnotatedBases(read, duplexIndels);

        expectedBases.get(7).setQual((byte) 0);
        expectedBases.get(8).setQual((byte) 0);

        expectedBases.get(5).deleteBase();
        expectedBases.get(6).deleteBase();

        expectedBases.get(5).setQual((byte) 0);
        expectedBases.get(6).setQual((byte) 0);

        assertTrue(readModified);
        assertEquals(expectedBases, annotedBases);
    }

    @Test
    public void testStripDuplexIndelsForwardRead()
    {
        int alignmentStart = 25;

        String refBases = "A".repeat(alignmentStart - 1 + 5) + "CT" + "A".repeat(1000);
        int refStart = 1;

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

        stripDuplexIndels(read, refBases.getBytes(), refStart);

        int alignmentScoreDiff = BWA_GAP_OPEN_PENALTY + 2 * BWA_GAP_EXTEND_PENALTY;

        int expectedMapq = min(MAX_MAPQ, round(mapq + ALIGNMENT_SCORED_DIFF_TO_MAPQ_DIFF * alignmentScoreDiff));
        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(5) + "CT" + "A".repeat(5);
        String expectedCigar = "12M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(5) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);

        SAMRecord expectedRead =
                createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, expectedReadStr, expectedCigar, false, false, null);

        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(expectedMapq);
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
        int refStart = 1;

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

        stripDuplexIndels(read, refBases.getBytes(), refStart);

        int alignmentScoreDiff = BWA_GAP_OPEN_PENALTY + 2 * BWA_GAP_EXTEND_PENALTY;

        int expectedMapq = min(MAX_MAPQ, round(mapq + ALIGNMENT_SCORED_DIFF_TO_MAPQ_DIFF * alignmentScoreDiff));
        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(5) + "CT" + "A".repeat(5);
        String expectedCigar = "12M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(5) + ZERO_QUAL.repeat(2) + NON_ZERO_QUAL.repeat(5);

        SAMRecord expectedRead =
                createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, expectedReadStr, expectedCigar, true, false, null);

        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(expectedMapq);
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
        int refStart = 1;

        fillQualZeroMismatchesWithRef(read, refBases.getBytes(), refStart);

        int alignmentScoreDiff = BWA_MISMATCH_PENALTY + BWA_MATCH_SCORE;

        int expectedMapq = min(MAX_MAPQ, round(mapq + ALIGNMENT_SCORED_DIFF_TO_MAPQ_DIFF * alignmentScoreDiff));
        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(10) + "G" + "A".repeat(9);
        String expectedCigar = "20M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(10) + phredToFastq(1) + NON_ZERO_QUAL.repeat(9);

        SAMRecord expectedRead =
                createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, expectedReadStr, expectedCigar, false, false, null);
        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(expectedMapq);
        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, expectedNm);
        expectedRead.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, expectedAlignmentScore);

        assertEquals(expectedRead, read);
    }
}
