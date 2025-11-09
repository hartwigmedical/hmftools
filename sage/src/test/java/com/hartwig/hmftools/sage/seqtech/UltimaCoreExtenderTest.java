package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsToStr;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.INVALID_INDEX;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.MISSING_BASE;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.addPadding;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.alignReadBases;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.extendCore;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.extendUltimaCore;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.extendUltimaCoreNew;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.findFlankBoundary;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.populateAlignedBaseLookupMaps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
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
import com.hartwig.hmftools.common.bam.ReadCigarState;
import com.hartwig.hmftools.common.utils.IntPair;
import com.hartwig.hmftools.sage.common.ReadCigarInfo;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.AlignedBase;
import com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.UltimaCoreInfo;

import org.junit.Test;

import htsjdk.samtools.CigarElement;

public class UltimaCoreExtenderTest
{
    private static final int TEST_FLANK_SIZE = 5;

    @Test
    public void testCoreExtensionBasic()
    {
        // pos / index                10             20             30        40
        //                  0123456789012345     67890     12345678901234567890123456789
        String refBases =  "AACCGGTTAACCGGTT" + "TACAT" + "TTAACCGGTTAACCGGTT";
        String readBases = refBases.substring(0, 17) + "G" + refBases.substring(18);

        RefSequence refSequence = new RefSequence(100, refBases.getBytes());

        int readAlignmentStart = 100;

        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(39, M));

        ReadCigarInfo readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, cigarElements, 106, 130, 116, 120,
                6, 30);

        ReadCigarInfo newReadCigarInfo = extendUltimaCoreNew(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(3, newReadCigarInfo.FlankIndexStart);
        assertEquals(103, newReadCigarInfo.FlankPositionStart);
        assertEquals(113, newReadCigarInfo.CorePositionStart);
        assertEquals(123, newReadCigarInfo.CorePositionEnd);
        assertEquals(33, newReadCigarInfo.FlankIndexEnd);
        assertEquals(133, newReadCigarInfo.FlankPositionEnd);
        assertEquals("31M", cigarElementsToStr(newReadCigarInfo.Cigar));

        // extend further due to non-matching bases from SNVs
        // pos / index         10             20             30        40
        //           0123456789012345     67890     12345678901234567890123456789
        refBases =  "AACCGGTTAACCAGGT" + "TACAT" + "TTTACCGGTTAACCGGTT";
        readBases = "AACCGGTTAACCGTTT" + "TAGAT" + "TTAGCCGGTTAACCGGTT"; // SNV at start and longer HP, then a shorter HP and another SNV

        refSequence = new RefSequence(100, refBases.getBytes());

        newReadCigarInfo = extendUltimaCoreNew(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(1, newReadCigarInfo.FlankIndexStart);
        assertEquals(101, newReadCigarInfo.FlankPositionStart);
        assertEquals(111, newReadCigarInfo.CorePositionStart);
        assertEquals(125, newReadCigarInfo.CorePositionEnd);
        assertEquals(35, newReadCigarInfo.FlankIndexEnd);
        assertEquals(135, newReadCigarInfo.FlankPositionEnd);
        assertEquals("35M", cigarElementsToStr(newReadCigarInfo.Cigar));
    }

    @Test
    public void testCoreExtensionIndels()
    {
        // pos / index                10             20             30        40
        //                  0123456789012345     67890     12345678901234567890123456789
        String refBases =  "AACCGGTTAACCGGTT" + "TACAT" + "TTAACCGGTTAACCGGTT";


        // pos / index                10  I          20         D        30        40
        //                  01234567890123456     78901     23     45678901234567890123456789
        String readBases = "AACCGGTTAACCGGGTT" + "TAGAT" + "TT" + "ACCGGTTAACCGGTT"; // has inserted G before core and deleted A after core

        RefSequence refSequence = new RefSequence(100, refBases.getBytes());

        int readAlignmentStart = 100;

        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(14, M),
                new CigarElement(1, I),
                new CigarElement(9, M),
                new CigarElement(1, D),
                new CigarElement(15, M));

        ReadCigarInfo readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, cigarElements, 106, 130, 116, 120,
                6, 30);

        ReadCigarInfo newReadCigarInfo = extendUltimaCoreNew(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(3, newReadCigarInfo.FlankIndexStart);
        assertEquals(103, newReadCigarInfo.FlankPositionStart);
        assertEquals(113, newReadCigarInfo.CorePositionStart);
        assertEquals(124, newReadCigarInfo.CorePositionEnd);
        assertEquals(34, newReadCigarInfo.FlankIndexEnd);
        assertEquals(134, newReadCigarInfo.FlankPositionEnd);
        assertEquals("11M1I9M1D11M", cigarElementsToStr(newReadCigarInfo.Cigar));
    }

    @Test
    public void testCoreExtensionSoftClips()
    {
        // pos / index                10             20             30        40
        //                  0123456789012345     67890     12345678901234567890123456789
        String refBases =  "AACCGGTTAACCGGTT" + "TACAT" + "TTAACCGGTTAACCGGTT";
        String readBases = refBases.substring(0, 17) + "G" + refBases.substring(18);

        RefSequence refSequence = new RefSequence(100, refBases.getBytes());

        int readAlignmentStart = 108;

        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(8, S),
                new CigarElement(31, M));

        ReadCigarInfo readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, cigarElements, 106, 130, 116, 120,
                6, 30);

        ReadCigarInfo newReadCigarInfo = extendUltimaCoreNew(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNull(newReadCigarInfo);

        // and now on the right in the core
        readAlignmentStart = 100;

        cigarElements = Lists.newArrayList(
                new CigarElement(22, M),
                new CigarElement(17, S));

        newReadCigarInfo = extendUltimaCoreNew(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNull(newReadCigarInfo);

        // invalid in the flanks
        cigarElements = Lists.newArrayList(
                new CigarElement(28, M),
                new CigarElement(11, S));

        newReadCigarInfo = extendUltimaCoreNew(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNull(newReadCigarInfo);

        // accept truncated flank in append mode
        cigarElements = Lists.newArrayList(
                new CigarElement(28, M),
                new CigarElement(11, S));

        newReadCigarInfo = extendUltimaCoreNew(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, true);

        assertNotNull(newReadCigarInfo);
        assertEquals(3, newReadCigarInfo.FlankIndexStart);
        assertEquals(103, newReadCigarInfo.FlankPositionStart);
        assertEquals(113, newReadCigarInfo.CorePositionStart);
        assertEquals(123, newReadCigarInfo.CorePositionEnd);
        assertEquals(27, newReadCigarInfo.FlankIndexEnd);
        assertEquals(127, newReadCigarInfo.FlankPositionEnd);
        assertEquals("25M", cigarElementsToStr(newReadCigarInfo.Cigar));
    }

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
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "A".repeat(1) + "G" + "A".repeat(9)
                + "T".repeat(2 * TEST_FLANK_SIZE);

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

        ReadCigarInfo readCigarInfo = new ReadCigarInfo(
                -1, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

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
        String refBases = "T".repeat(2 * TEST_FLANK_SIZE) + "GAA" + "A".repeat(1) + "G" + "A".repeat(9)
                + "AAG" + "T".repeat(2 * TEST_FLANK_SIZE);

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

        ReadCigarInfo readCigarInfo = new ReadCigarInfo(
                4, null, -1, -1, corePositionStart, corePositionEnd, -1, -1);

        UltimaCoreInfo ultimaCoreInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, TEST_FLANK_SIZE, false);

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
