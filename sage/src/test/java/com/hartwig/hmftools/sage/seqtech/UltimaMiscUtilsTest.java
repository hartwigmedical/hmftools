package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsToStr;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.sage.common.TestUtils.setUltimaSequencing;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.extendUltimaCore;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.isAdjacentToLongHomopolymer;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.isMsiIndelOfType;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.ReadCigarInfo;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class UltimaMiscUtilsTest
{
    public UltimaMiscUtilsTest()
    {
        setUltimaSequencing();
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }

    @Test
    public void testIsMsiIndelTypes()
    {
        // not an indel
        SimpleVariant variant = new SimpleVariant(CHR_1, 100, "A", "C");
        assertFalse(isMsiIndelOfType(variant, null));

        // no repeat units
        variant = new SimpleVariant(CHR_1, 100, "A", "AA");
        List<String> units = Lists.newArrayList();
        assertFalse(isMsiIndelOfType(variant, units));

        // single unit but not MSI
        variant = new SimpleVariant(CHR_1, 100, "A", "AAB");
        units = Lists.newArrayList("A");
        assertFalse(isMsiIndelOfType(variant, units));

        // MSI repeat
        variant = new SimpleVariant(CHR_1, 100, "A", "AABABAB");
        units = Lists.newArrayList("AB");
        assertTrue(isMsiIndelOfType(variant, units));
    }

    @Test
    public void testIsAdjacentToLongHomopolymers()
    {
        int longLength = 10;
        String readBases = "A".repeat(2*longLength - 1);
        SAMRecord read = createRead(readBases);

        // not long enough
        assertFalse(isAdjacentToLongHomopolymer(read, longLength - 1, longLength));

        readBases = "A".repeat(longLength) + "C" + "T".repeat(longLength - 1);
        read = createRead(readBases);

        // on left
        assertTrue(isAdjacentToLongHomopolymer(read, longLength, longLength));

        readBases = "T".repeat(longLength - 1) + "C" + "A".repeat(longLength);
        read = createRead(readBases);

        // on right
        assertTrue(isAdjacentToLongHomopolymer(read, longLength - 1, longLength));

        readBases = "C" + "A".repeat(2*longLength - 1);
        read = createRead(readBases);

        // not on left
        assertFalse(isAdjacentToLongHomopolymer(read, longLength, longLength));

        readBases = "A".repeat(2*longLength - 1) + "C";
        read = createRead(readBases);

        // not on right
        assertFalse(isAdjacentToLongHomopolymer(read, longLength - 1, longLength));

        readBases = "N".repeat(2*longLength);
        read = createRead(readBases);

        // invalid bases
        assertFalse(isAdjacentToLongHomopolymer(read, longLength, longLength));
    }

    private static SAMRecord createRead(final String readBases)
    {
        String cigar = format("%dM", readBases.length());
        return buildSamRecord(100, cigar, readBases);
    }

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

        ReadCigarInfo newReadCigarInfo = extendUltimaCore(
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

        newReadCigarInfo = extendUltimaCore(
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

        ReadCigarInfo newReadCigarInfo = extendUltimaCore(
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

        ReadCigarInfo newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNull(newReadCigarInfo);

        // and now on the right in the core
        readAlignmentStart = 100;

        cigarElements = Lists.newArrayList(
                new CigarElement(22, M),
                new CigarElement(17, S));

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNull(newReadCigarInfo);

        // invalid in the flanks
        cigarElements = Lists.newArrayList(
                new CigarElement(28, M),
                new CigarElement(11, S));

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNull(newReadCigarInfo);

        // accept truncated flank in append mode
        cigarElements = Lists.newArrayList(
                new CigarElement(28, M),
                new CigarElement(11, S));

        newReadCigarInfo = extendUltimaCore(
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
    public void testBuildContextUltimaNoCoreExtension()
    {
        String middleRefBases = "CTCGGCCTCCCGGAGGTGCCGG";
        String middleReadBases = "CTCGGCCTCCCTGAGGTGCCGG";
        int variantMiddleIndex = 11;
        int refPaddingSize = 20;
        int readPaddingSize = 20;
        String refPadding = "A".repeat(refPaddingSize);
        String readPadding = "A".repeat(readPaddingSize);
        String refBases = refPadding + middleRefBases + refPadding;
        String readBases = readPadding + middleReadBases + readPadding;

        SimpleVariant variant = new SimpleVariant(CHR_1, 1 + refPaddingSize + variantMiddleIndex, "G", "T");

        int alignmentStart = 1 + refPaddingSize - readPaddingSize;
        String cigar = readBases.length() + "M";
        SAMRecord read = createSamRecordUnpaired(
                "READ_001", CHR_1, alignmentStart, readBases, cigar, false, false, null);

        int varIndexInRead = readPaddingSize + variantMiddleIndex;

        RefSequence refSequence = new RefSequence(1, refBases.getBytes());

        // SageConfig config = new SageConfig(false, SequencingType.ULTIMA);
        VariantReadContextBuilder readContextBuilder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = readContextBuilder.buildContext(variant, read, varIndexInRead, refSequence);
        int expectedCoreLength = 7;

        assertNotNull(readContext);
        assertEquals(DEFAULT_FLANK_LENGTH, readContext.CoreIndexStart);
        assertEquals(DEFAULT_FLANK_LENGTH + expectedCoreLength - 1, readContext.CoreIndexEnd);
        assertEquals(variant.Position - 4, readContext.CorePositionStart);
        assertEquals(variant.Position + 2, readContext.CorePositionEnd);
        assertEquals("TCCCGGA", readContext.refBases());
        assertEquals("TCCCTGA", readContext.coreStr());
        assertEquals(
                readBases.substring(varIndexInRead - 4 - DEFAULT_FLANK_LENGTH, varIndexInRead + 3 + DEFAULT_FLANK_LENGTH),
                readContext.readBases());
    }

    @Test
    public void testBuildContextUltimaWithExtension()
    {
        String middleRefBases = "CTTTTTTTTTTAATTGCAA";
        String middleReadBases = "CTTTTTTTTTAATTGGCAA";
        int variantMiddleIndex = 14;
        int refPaddingSize = 20;
        int readPaddingSize = 20;
        String refBases = "A".repeat(refPaddingSize) + middleRefBases + "T".repeat(refPaddingSize);
        String readBases = "A".repeat(readPaddingSize) + middleReadBases + "T".repeat(readPaddingSize);

        SimpleVariant variant = new SimpleVariant(CHR_1, 1 + refPaddingSize + variantMiddleIndex, "T", "G");

        int alignmentStart = 1 + refPaddingSize - readPaddingSize;
        String cigar = readBases.length() + "M";
        SAMRecord read = createSamRecordUnpaired(
                "READ_001", CHR_1, alignmentStart, readBases, cigar, false, false, null);

        int varIndexInRead = readPaddingSize + variantMiddleIndex;

        RefSequence refSequence = new RefSequence(1, refBases.getBytes());

        VariantReadContextBuilder readContextBuilder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = readContextBuilder.buildContext(variant, read, varIndexInRead, refSequence);
        int expectedCoreLength = 17;

        assertNotNull(readContext);
        assertEquals(DEFAULT_FLANK_LENGTH, readContext.CoreIndexStart);
        assertEquals(DEFAULT_FLANK_LENGTH + expectedCoreLength - 1, readContext.CoreIndexEnd);
        assertEquals(variant.Position - 14, readContext.CorePositionStart);
        assertEquals(variant.Position + 2, readContext.CorePositionEnd);
        assertEquals("CTTTTTTTTTTAATTGC", readContext.refBases());
        assertEquals("CTTTTTTTTTAATTGGC", readContext.coreStr());
        assertEquals(
                readBases.substring(varIndexInRead - 14 - DEFAULT_FLANK_LENGTH, varIndexInRead + 3 + DEFAULT_FLANK_LENGTH),
                readContext.readBases());
    }
}
