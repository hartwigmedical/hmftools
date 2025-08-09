package com.hartwig.hmftools.sage.ultima;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.INVALID_INDEX;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.MISSING_BASE;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.addPadding;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.alignReadBases;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.extendCore;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.extendUltimaCore;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.findFlankBoundary;
import static com.hartwig.hmftools.sage.common.UltimaCoreExtender.populateAlignedBaseLookupMaps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.IntPair;
import com.hartwig.hmftools.sage.common.ReadCigarInfo;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.UltimaCoreExtender.AlignedBase;
import com.hartwig.hmftools.sage.common.UltimaCoreExtender.UltimaCoreInfo;

import org.junit.Test;

import htsjdk.samtools.CigarElement;

public class UltimaCoreExtenderTest
{
    private static final int TEST_FLANK_SIZE = 5;

    @Test
    public void testAlignReadBasesWithIndels()
    {
        String refBases = "AAAAGAA";
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        String readBases = "AAGAAAA";
        int readAlignmentStart = 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2, M),
                new CigarElement(1, I),
                new CigarElement(2, M),
                new CigarElement(1, D),
                new CigarElement(2, M));

        List<AlignedBase> actualAlignment = alignReadBases(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements);
        List<AlignedBase> expectedAlignment = Lists.newArrayList(
                new AlignedBase(1, 0, (byte) 'A', (byte) 'A', M),
                new AlignedBase(2, 1, (byte) 'A', (byte) 'A', M),
                new AlignedBase(2, 2, MISSING_BASE, (byte) 'G', I),
                new AlignedBase(3, 3, (byte) 'A', (byte) 'A', M),
                new AlignedBase(4, 4, (byte) 'A', (byte) 'A', M),
                new AlignedBase(5, 4, (byte) 'G', MISSING_BASE, D),
                new AlignedBase(6, 5, (byte) 'A', (byte) 'A', M),
                new AlignedBase(7, 6, (byte) 'A', (byte) 'A', M));

        assertEquals(expectedAlignment, actualAlignment);
    }

    @Test
    public void testAlignReadBasesLeftHardClip()
    {
        String refBases = "A";
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        String readBases = "A";
        int readAlignmentStart = 1;
        List<CigarElement> cigarElements = Lists.newArrayList(new CigarElement(1, H), new CigarElement(1, M));

        List<AlignedBase> actualAlignment = alignReadBases(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements);
        List<AlignedBase> expectedAlignment = Lists.newArrayList(new AlignedBase(1, 0, (byte) 'A', (byte) 'A', M));

        assertEquals(expectedAlignment, actualAlignment);
    }

    @Test
    public void testAlignReadBasesIgnoreRightSoftClip()
    {
        String refBases = "AG";
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        String readBases = "AT";
        int readAlignmentStart = 1;
        List<CigarElement> cigarElements = Lists.newArrayList(new CigarElement(1, M), new CigarElement(1, S));

        List<AlignedBase> actualAlignment = alignReadBases(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements);
        List<AlignedBase> expectedAlignment = Lists.newArrayList(new AlignedBase(1, 0, (byte) 'A', (byte) 'A', M));

        assertEquals(expectedAlignment, actualAlignment);
    }

    @Test
    public void testAlignReadBasesCaptureLeftSoftClip()
    {
        String refBases = "GA";
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        String readBases = "TA";
        int readAlignmentStart = 2;
        List<CigarElement> cigarElements = Lists.newArrayList(new CigarElement(1, S), new CigarElement(1, M));

        List<AlignedBase> actualAlignment = alignReadBases(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements);
        List<AlignedBase> expectedAlignment = Lists.newArrayList(
                new AlignedBase(1, 0, (byte) 'G', (byte) 'T', S),
                new AlignedBase(2, 1, (byte) 'A', (byte) 'A', M));

        assertEquals(expectedAlignment, actualAlignment);
    }

    @Test
    public void testExtendCoreBasesAtLeftAllMissing()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, MISSING_BASE, (byte) 'A', I),
                new AlignedBase(2, 1, (byte) 'A', (byte) 'A', M));
        Map<Integer, Integer> lookupFromRefPos = Maps.newHashMap();
        Map<Integer, Integer> lookupFromReadIndex = Maps.newHashMap();
        populateAlignedBaseLookupMaps(alignedBases, lookupFromRefPos, lookupFromReadIndex);

        IntPair extendCoreIndices = extendCore(alignedBases, 0, 1, lookupFromRefPos, lookupFromReadIndex);

        assertNull(extendCoreIndices);
    }

    @Test
    public void testExtendCoreBasesAtRightAllMissing()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, (byte) 'A', (byte) 'A', M),
                new AlignedBase(1, 1, MISSING_BASE, (byte) 'A', I));
        Map<Integer, Integer> lookupFromRefPos = Maps.newHashMap();
        Map<Integer, Integer> lookupFromReadIndex = Maps.newHashMap();
        populateAlignedBaseLookupMaps(alignedBases, lookupFromRefPos, lookupFromReadIndex);

        IntPair extendCoreIndices = extendCore(alignedBases, 0, 1, lookupFromRefPos, lookupFromReadIndex);

        assertNull(extendCoreIndices);
    }

    @Test
    public void testExtendCoreExtendToCoverReadHomopolyerAtEndOfCore()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, (byte) 'G', (byte) 'G', M),
                new AlignedBase(2, 1, (byte) 'A', (byte) 'A', M),
                new AlignedBase(3, 2, (byte) 'T', (byte) 'A', M),
                new AlignedBase(4, 3, (byte) 'A', (byte) 'A', M),
                new AlignedBase(5, 4, (byte) 'G', (byte) 'G', M)
        );
        Map<Integer, Integer> lookupFromRefPos = Maps.newHashMap();
        Map<Integer, Integer> lookupFromReadIndex = Maps.newHashMap();
        populateAlignedBaseLookupMaps(alignedBases, lookupFromRefPos, lookupFromReadIndex);

        IntPair extendCoreIndices = extendCore(alignedBases, 1, 2, lookupFromRefPos, lookupFromReadIndex);
        int coreStart = extendCoreIndices.getLeft();
        int coreEnd = extendCoreIndices.getRight();

        assertEquals(1, coreStart);
        assertEquals(3, coreEnd);
    }

    @Test
    public void testExtendCoreInitCoreIndicesOutOfRange()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, DNA_BASE_BYTES[0], DNA_BASE_BYTES[0], M),
                new AlignedBase(2, 1, DNA_BASE_BYTES[0], DNA_BASE_BYTES[0], M),
                new AlignedBase(3, 2, DNA_BASE_BYTES[0], DNA_BASE_BYTES[0], M));

        IntPair extendCoreIndices = extendCore(alignedBases, -1, 2, null, null);
        assertNull(extendCoreIndices);

        extendCoreIndices = extendCore(alignedBases, 1, 3, null, null);
        assertNull(extendCoreIndices);
    }

    @Test
    public void testAddPaddingNoSpaceForPadding()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, DNA_BASE_BYTES[0], DNA_BASE_BYTES[0], M));

        int coreStart = addPadding(alignedBases, 0, true);
        assertEquals(INVALID_INDEX, coreStart);

        int coreEnd = addPadding(alignedBases, 1, false);
        assertEquals(INVALID_INDEX, coreEnd);
    }

    @Test
    public void testAddPaddingPaddingAtNextBase()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, (byte) 'A', (byte) 'A', M),
                new AlignedBase(2, 1, (byte) 'A', (byte) 'A', M));

        int coreStart = addPadding(alignedBases, 1, true);
        assertEquals(0, coreStart);

        int coreEnd = addPadding(alignedBases, 0, false);
        assertEquals(1, coreEnd);
    }

    @Test
    public void testFindFlankBoundaryNotEnoughBases()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, MISSING_BASE, (byte) 'A', I),
                new AlignedBase(2, 1, (byte) 'A', (byte) 'A', M),
                new AlignedBase(2, 2, MISSING_BASE, (byte) 'A', I));

        int flankStart = findFlankBoundary(alignedBases, 1, true, 2, false);
        assertEquals(INVALID_INDEX, flankStart);

        int flankEnd = findFlankBoundary(alignedBases, 1, false, 2, false);
        assertEquals(INVALID_INDEX, flankEnd);
    }

    @Test
    public void testFindFlankBoundaryPushAlignmentOut()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, (byte) 'A', (byte) 'A', M),
                new AlignedBase(1, 1, MISSING_BASE, (byte) 'A', I),
                new AlignedBase(2, 2, (byte) 'A', (byte) 'A', M),
                new AlignedBase(2, 3, MISSING_BASE, (byte) 'A', I),
                new AlignedBase(3, 4, (byte) 'A', (byte) 'A', M));

        int flankStart = findFlankBoundary(alignedBases, 2, true, 1, false);
        assertEquals(0, flankStart);

        int flankEnd = findFlankBoundary(alignedBases, 2, false, 1, false);
        assertEquals(4, flankEnd);
    }

    @Test
    public void testFindFlankBoundaryPushAlignmentOutExhaustsBases()
    {
        List<AlignedBase> alignedBases = Lists.newArrayList(
                new AlignedBase(1, 0, MISSING_BASE, (byte) 'A', I),
                new AlignedBase(2, 1, (byte) 'A', (byte) 'A', M),
                new AlignedBase(2, 2, MISSING_BASE, (byte) 'A', I));

        int flankStart = findFlankBoundary(alignedBases, 1, true, 1, false);
        assertEquals(INVALID_INDEX, flankStart);

        int flankEnd = findFlankBoundary(alignedBases, 1, false, 1, false);
        assertEquals(INVALID_INDEX, flankEnd);
    }

    @Test
    public void testExtendUltimaCoreNoCoreExtensionWithInsert()
    {
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(9) + "T".repeat(2 * TEST_FLANK_SIZE);
        String readBases = "C".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(1) + "G" + "A".repeat(8) + "T".repeat(2 * TEST_FLANK_SIZE);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        int readAlignmentStart = 2 * TEST_FLANK_SIZE + 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2 * TEST_FLANK_SIZE, S),
                new CigarElement(1, M),
                new CigarElement(1, I),
                new CigarElement(8 + 2 * TEST_FLANK_SIZE, M));
        int corePositionStart = 2 * TEST_FLANK_SIZE + 1;
        int corePositionEnd = 2 * TEST_FLANK_SIZE + 9;
        ReadCigarInfo readCigarInfo = new ReadCigarInfo(-1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

        assertEquals(2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreStart);
        assertEquals(readBases.length() - 1 - 2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreEnd);
        assertTrue(ultimaCoreInfo.CigarInfo.isValid());
        assertEquals(
                Lists.newArrayList(
                        new CigarElement(TEST_FLANK_SIZE, S),
                        new CigarElement(1, M),
                        new CigarElement(1, I),
                        new CigarElement(8 + TEST_FLANK_SIZE, M)),
                ultimaCoreInfo.CigarInfo.Cigar);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.ReadAlignmentStart);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.FlankPositionStart);
        assertEquals(corePositionEnd + TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankPositionEnd);
        assertEquals(corePositionStart, ultimaCoreInfo.CigarInfo.CorePositionStart);
        assertEquals(corePositionEnd, ultimaCoreInfo.CigarInfo.CorePositionEnd);
        assertEquals(TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexStart);
        assertEquals(readBases.length() - 1 - TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexEnd);
    }

    @Test
    public void testExtendUltimaCoreNoCoreExtensionWithDel()
    {
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(1) + "G" + "A".repeat(9) + "T".repeat(2 * TEST_FLANK_SIZE);
        String readBases = "C".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(10) + "T".repeat(2 * TEST_FLANK_SIZE);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        int readAlignmentStart = 2 * TEST_FLANK_SIZE + 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2 * TEST_FLANK_SIZE, S),
                new CigarElement(1, M),
                new CigarElement(1, D),
                new CigarElement(9 + 2 * TEST_FLANK_SIZE, M));
        int corePositionStart = 2 * TEST_FLANK_SIZE + 1;
        int corePositionEnd = 2 * TEST_FLANK_SIZE + 11;
        ReadCigarInfo readCigarInfo = new ReadCigarInfo(-1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

        assertEquals(2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreStart);
        assertEquals(readBases.length() - 1 - 2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreEnd);
        assertTrue(ultimaCoreInfo.CigarInfo.isValid());
        assertEquals(
                Lists.newArrayList(
                        new CigarElement(TEST_FLANK_SIZE, S),
                        new CigarElement(1, M),
                        new CigarElement(1, D),
                        new CigarElement(9 + TEST_FLANK_SIZE, M)),
                ultimaCoreInfo.CigarInfo.Cigar);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.ReadAlignmentStart);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.FlankPositionStart);
        assertEquals(corePositionEnd + TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankPositionEnd);
        assertEquals(corePositionStart, ultimaCoreInfo.CigarInfo.CorePositionStart);
        assertEquals(corePositionEnd, ultimaCoreInfo.CigarInfo.CorePositionEnd);
        assertEquals(TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexStart);
        assertEquals(readBases.length() - 1 - TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexEnd);
    }

    @Test
    public void testExtendUltimaCore()
    {
        String refBases =
                "T".repeat(2 * TEST_FLANK_SIZE) + "GAA" + "A".repeat(1) + "G" + "A".repeat(9) + "AAG" + "T".repeat(2 * TEST_FLANK_SIZE);
        String readBases = "C".repeat(2 * TEST_FLANK_SIZE) + "GAA" + "A".repeat(10) + "AAG" + "T".repeat(2 * TEST_FLANK_SIZE);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        int readAlignmentStart = 2 * TEST_FLANK_SIZE + 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2 * TEST_FLANK_SIZE, S),
                new CigarElement(4, M),
                new CigarElement(1, D),
                new CigarElement(12 + 2 * TEST_FLANK_SIZE, M));
        int corePositionStart = 2 * TEST_FLANK_SIZE + 4;
        int corePositionEnd = 2 * TEST_FLANK_SIZE + 4 + 10;
        ReadCigarInfo readCigarInfo = new ReadCigarInfo(-1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

        assertEquals(2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreStart);
        assertEquals(readBases.length() - 1 - 2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreEnd);
        assertTrue(ultimaCoreInfo.CigarInfo.isValid());
        assertEquals(
                Lists.newArrayList(
                        new CigarElement(TEST_FLANK_SIZE, S),
                        new CigarElement(4, M),
                        new CigarElement(1, D),
                        new CigarElement(12 + TEST_FLANK_SIZE, M)),
                ultimaCoreInfo.CigarInfo.Cigar);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.ReadAlignmentStart);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.FlankPositionStart);
        assertEquals(corePositionEnd + 3 + TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankPositionEnd);
        assertEquals(corePositionStart - 3, ultimaCoreInfo.CigarInfo.CorePositionStart);
        assertEquals(corePositionEnd + 3, ultimaCoreInfo.CigarInfo.CorePositionEnd);
        assertEquals(TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexStart);
        assertEquals(readBases.length() - 1 - TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexEnd);
    }

    @Test
    public void testExtendUltimaCoreExtendedCoreMustBeginOnMatchInserts()
    {
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "AGA" + "A".repeat(10) + "T".repeat(2 * TEST_FLANK_SIZE);
        String readBases = "C".repeat(2 * TEST_FLANK_SIZE) + "AGAA" + "A".repeat(10) + "T".repeat(2 * TEST_FLANK_SIZE);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        int readAlignmentStart = 2 * TEST_FLANK_SIZE + 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2 * TEST_FLANK_SIZE, S),
                new CigarElement(2, M),
                new CigarElement(1, I),
                new CigarElement(11 + 2 * TEST_FLANK_SIZE, M));
        int corePositionStart = 2 * TEST_FLANK_SIZE + 4;
        int corePositionEnd = 2 * TEST_FLANK_SIZE + 13;
        ReadCigarInfo readCigarInfo = new ReadCigarInfo(-1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

        assertEquals(2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreStart);
        assertEquals(readBases.length() - 1 - 2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreEnd);
        assertTrue(ultimaCoreInfo.CigarInfo.isValid());
        assertEquals(
                Lists.newArrayList(
                        new CigarElement(TEST_FLANK_SIZE, S),
                        new CigarElement(2, M),
                        new CigarElement(1, I),
                        new CigarElement(11 + TEST_FLANK_SIZE, M)),
                ultimaCoreInfo.CigarInfo.Cigar);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.ReadAlignmentStart);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.FlankPositionStart);
        assertEquals(corePositionEnd + TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankPositionEnd);
        assertEquals(corePositionStart - 3, ultimaCoreInfo.CigarInfo.CorePositionStart);
        assertEquals(corePositionEnd, ultimaCoreInfo.CigarInfo.CorePositionEnd);
        assertEquals(TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexStart);
        assertEquals(readBases.length() - 1 - TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexEnd);
    }

    @Test
    public void testExtendUltimaCoreExtendedCoreMustEndOnMatchInserts()
    {
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(9) + "GGTA" + "T".repeat(2 * TEST_FLANK_SIZE);
        String readBases = "C".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(9) + "GGGTA" + "T".repeat(2 * TEST_FLANK_SIZE);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        int readAlignmentStart = 2 * TEST_FLANK_SIZE + 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2 * TEST_FLANK_SIZE, S),
                new CigarElement(9, M),
                new CigarElement(1, I),
                new CigarElement(4 + 2 * TEST_FLANK_SIZE, M));
        int corePositionStart = 2 * TEST_FLANK_SIZE + 1;
        int corePositionEnd = 2 * TEST_FLANK_SIZE + 10;
        ReadCigarInfo readCigarInfo = new ReadCigarInfo(-1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

        assertEquals(2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreStart);
        assertEquals(readBases.length() - 2 - 2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreEnd);
        assertTrue(ultimaCoreInfo.CigarInfo.isValid());
        assertEquals(
                Lists.newArrayList(
                        new CigarElement(TEST_FLANK_SIZE, S),
                        new CigarElement(9, M),
                        new CigarElement(1, I),
                        new CigarElement(3 + TEST_FLANK_SIZE, M)),
                ultimaCoreInfo.CigarInfo.Cigar);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.ReadAlignmentStart);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.FlankPositionStart);
        assertEquals(corePositionEnd + 2 + TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankPositionEnd);
        assertEquals(corePositionStart, ultimaCoreInfo.CigarInfo.CorePositionStart);
        assertEquals(corePositionEnd + 2, ultimaCoreInfo.CigarInfo.CorePositionEnd);
        assertEquals(TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexStart);
        assertEquals(readBases.length() - 2 - TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexEnd);
    }

    @Test
    public void testExtendUltimaCoreExtendedCoreMustBeginOnMatchDels()
    {
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "AGAA" + "A".repeat(10) + "T".repeat(2 * TEST_FLANK_SIZE);
        String readBases = "C".repeat(2 * TEST_FLANK_SIZE) + "AGA" + "A".repeat(10) + "T".repeat(2 * TEST_FLANK_SIZE);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        int readAlignmentStart = 2 * TEST_FLANK_SIZE + 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2 * TEST_FLANK_SIZE, S),
                new CigarElement(2, M),
                new CigarElement(1, D),
                new CigarElement(11 + 2 * TEST_FLANK_SIZE, M));
        int corePositionStart = 2 * TEST_FLANK_SIZE + 5;
        int corePositionEnd = 2 * TEST_FLANK_SIZE + 14;
        ReadCigarInfo readCigarInfo = new ReadCigarInfo(-1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

        assertEquals(2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreStart);
        assertEquals(readBases.length() - 1 - 2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreEnd);
        assertTrue(ultimaCoreInfo.CigarInfo.isValid());
        assertEquals(
                Lists.newArrayList(
                        new CigarElement(TEST_FLANK_SIZE, S),
                        new CigarElement(2, M),
                        new CigarElement(1, D),
                        new CigarElement(11 + TEST_FLANK_SIZE, M)),
                ultimaCoreInfo.CigarInfo.Cigar);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.ReadAlignmentStart);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.FlankPositionStart);
        assertEquals(corePositionEnd + TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankPositionEnd);
        assertEquals(corePositionStart - 4, ultimaCoreInfo.CigarInfo.CorePositionStart);
        assertEquals(corePositionEnd, ultimaCoreInfo.CigarInfo.CorePositionEnd);
        assertEquals(TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexStart);
        assertEquals(readBases.length() - 1 - TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexEnd);
    }

    @Test
    public void testExtendUltimaCoreExtendedCoreMustEndOnMatchDels()
    {
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(9) + "GGGTA" + "T".repeat(2 * TEST_FLANK_SIZE);
        String readBases = "C".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(9) + "GGTA" + "T".repeat(2 * TEST_FLANK_SIZE);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        int readAlignmentStart = 2 * TEST_FLANK_SIZE + 1;
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(2 * TEST_FLANK_SIZE, S),
                new CigarElement(9, M),
                new CigarElement(1, D),
                new CigarElement(4 + 2 * TEST_FLANK_SIZE, M));
        int corePositionStart = 2 * TEST_FLANK_SIZE + 1;
        int corePositionEnd = 2 * TEST_FLANK_SIZE + 10;
        ReadCigarInfo readCigarInfo = new ReadCigarInfo(-1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

        assertEquals(2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreStart);
        assertEquals(readBases.length() - 2 - 2 * TEST_FLANK_SIZE, ultimaCoreInfo.ReadCoreEnd);
        assertTrue(ultimaCoreInfo.CigarInfo.isValid());
        assertEquals(
                Lists.newArrayList(
                        new CigarElement(TEST_FLANK_SIZE, S),
                        new CigarElement(9, M),
                        new CigarElement(1, D),
                        new CigarElement(3 + TEST_FLANK_SIZE, M)),
                ultimaCoreInfo.CigarInfo.Cigar);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.ReadAlignmentStart);
        assertEquals(readAlignmentStart, ultimaCoreInfo.CigarInfo.FlankPositionStart);
        assertEquals(corePositionEnd + 3 + TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankPositionEnd);
        assertEquals(corePositionStart, ultimaCoreInfo.CigarInfo.CorePositionStart);
        assertEquals(corePositionEnd + 3, ultimaCoreInfo.CigarInfo.CorePositionEnd);
        assertEquals(TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexStart);
        assertEquals(readBases.length() - 2 - TEST_FLANK_SIZE, ultimaCoreInfo.CigarInfo.FlankIndexEnd);
    }

    @Test
    public void testExtendUltimaCoreLeftSoftClipOverRefGenomeBoundary()
    {
        String refBases = "GA";
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        String readBases = "TTA";
        int readAlignmentStart = 2;
        List<CigarElement> cigarElements = Lists.newArrayList(new CigarElement(2, S), new CigarElement(1, M));

        UltimaCoreInfo ultimaCoreInfo =
                extendUltimaCore(readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, null, INVALID_INDEX, false);

        assertNull(ultimaCoreInfo);
    }
}
