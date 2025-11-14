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
        // ref pos                    10                       20                      30        40
        //                  012345678901     I     2345     67890     12     3     45678901234567890123456789
        String refBases =  "AACCGGTTAACC" +  "" + "GGTT" + "TACAT" + "TT" + "A" + "ACCGGTTAACCGGTT";

        // read index                 10                       20                        30        40
        //                  012345678901     2     3456     78901     23     D     45678901234567890123456789
        String readBases = "AACCGGTTAACC" + "G" + "GGTT" + "TAGAT" + "TT" +  "" + "ACCGGTTAACCGGTT"; // has inserted G before core and deleted A after core

        RefSequence refSequence = new RefSequence(100, refBases.getBytes());

        int readAlignmentStart = 100;

        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(12, M),
                new CigarElement(1, I),
                new CigarElement(11, M),
                new CigarElement(1, D),
                new CigarElement(15, M));

        ReadCigarInfo readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, cigarElements, 107, 130, 116, 120,
                17, 30);

        ReadCigarInfo newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(4, newReadCigarInfo.FlankIndexStart);
        assertEquals(104, newReadCigarInfo.FlankPositionStart);
        assertEquals(113, newReadCigarInfo.CorePositionStart);
        assertEquals(124, newReadCigarInfo.CorePositionEnd);
        assertEquals(34, newReadCigarInfo.FlankIndexEnd);
        assertEquals(134, newReadCigarInfo.FlankPositionEnd);
        assertEquals("8M1I11M1D11M", cigarElementsToStr(newReadCigarInfo.Cigar));

        // a deleted base immediately before and after the core

        // pos                 10                       20                  30
        //           012345678901234     5     67     89012     3     456789012345678
        refBases =  "AACCGGTTAACCGGG" + "A" + "AA" + "TAGAT" + "A" + "ACCGGTTAACCGGTT";

        //                     10                        20                   30
        // index     012345678901234     D     56     78901           234567890123456
        readBases = "AACCGGTTAACCGGG" +  "" + "AA" + "TAGAT" +  "" + "ACCGGTTAACCGGTT";
        // flank          FFFFFFFFFF                                  FFFFFFFFFF

        refSequence = new RefSequence(100, refBases.getBytes());

        readAlignmentStart = 100;

        cigarElements = Lists.newArrayList(
                new CigarElement(15, M),
                new CigarElement(1, D),
                new CigarElement(7, M),
                new CigarElement(1, D),
                new CigarElement(15, M));

        List<CigarElement> readCigar = List.of(
                new CigarElement(8, M),
                new CigarElement(1, D),
                new CigarElement(7, M),
                new CigarElement(1, D),
                new CigarElement(10, M));

        readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, readCigar, 107, 133, 118, 122,
                7, 31);

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(7, newReadCigarInfo.FlankIndexStart);
        assertEquals(107, newReadCigarInfo.FlankPositionStart);
        assertEquals(118, newReadCigarInfo.CorePositionStart);
        assertEquals(122, newReadCigarInfo.CorePositionEnd);
        assertEquals(31, newReadCigarInfo.FlankIndexEnd);
        assertEquals(133, newReadCigarInfo.FlankPositionEnd);
        assertEquals("8M1D7M1D10M", cigarElementsToStr(newReadCigarInfo.Cigar));

        // an inserted base immediately before and after the core

        // pos                 10                  20        30
        //           012345678901234     56789     0123456789012345678
        refBases =  "AACCGGTTAACCAGG" + "ATCTA" + "TCCGGTTAACCGGTT";

        //                     10                               20        30
        // pos       012345678901234           56789            01234567890123456
        //                     10                  20                   30
        // index     012345678901234     5     67890      1     234567890123456
        readBases = "AACCGGTTAACCAGG" + "T" + "ATCTA" +  "T" + "TCCGGTTAACCGGTT";
        // flank          FFFFFFFFFF                            FFFFFFFFFF

        refSequence = new RefSequence(100, refBases.getBytes());

        readAlignmentStart = 100;

        cigarElements = Lists.newArrayList(
                new CigarElement(15, M),
                new CigarElement(1, I),
                new CigarElement(5, M),
                new CigarElement(1, I),
                new CigarElement(15, M));

        readCigar = List.of(
                new CigarElement(11, M),
                new CigarElement(1, I),
                new CigarElement(5, M),
                new CigarElement(1, I),
                new CigarElement(11, M));

        readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, readCigar, 105, 129, 115, 119,
                5, 31);

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertTrue(checkEqual(readCigarInfo, newReadCigarInfo));
    }

    @Test
    public void testCoreExtensionMiscExamples()
    {
        // test 1: HP starts on ref but switches to read then ref again

        // pos / index                10             20             30        40
        //                  0123456789012345     67890     12345678901234567890123456789
        String refBases =  "AACCGGTTAACCGGTA" + "ACTCA" + "ATAAACGGTTAACCGGTT";
        String readBases = "AACCGGTTAACCGTTT" + "ACTCA" + "TTTACCGGTTAACCGGTT";
        //                        FFFFFFFFFF               FFFFFFFFFF

        RefSequence refSequence = new RefSequence(100, refBases.getBytes());

        int readAlignmentStart = 100;

        List<CigarElement> cigarElements = Lists.newArrayList(new CigarElement(39, M));
        List<CigarElement> readCigarElements = Lists.newArrayList(new CigarElement(25, M));

        ReadCigarInfo readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, readCigarElements, 106, 130, 116, 120,
                6, 30);

        ReadCigarInfo newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(1, newReadCigarInfo.FlankIndexStart);
        assertEquals(101, newReadCigarInfo.FlankPositionStart);
        assertEquals(111, newReadCigarInfo.CorePositionStart);
        assertEquals(127, newReadCigarInfo.CorePositionEnd);
        assertEquals(37, newReadCigarInfo.FlankIndexEnd);
        assertEquals(137, newReadCigarInfo.FlankPositionEnd);
        assertEquals("37M", cigarElementsToStr(newReadCigarInfo.Cigar));

        // test 2: no HP, but next bases don't match and so require extension

        // pos / index         10             20             30        40
        //           0123456789012345     67890     12345678901234567890123456789
        refBases =  "AACCGGTTAACCGATC" + "ACTCA" + "CTAGACGGTTAACCGGTT";
        readBases = "AACCGGTTAACCGATT" + "GCTCG" + "TTAGACGGTTAACCGGTT";
        //                 FFFFFFFFFF               FFFFFFFFFF

        refSequence = new RefSequence(100, refBases.getBytes());

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(4, newReadCigarInfo.FlankIndexStart);
        assertEquals(104, newReadCigarInfo.FlankPositionStart);
        assertEquals(114, newReadCigarInfo.CorePositionStart);
        assertEquals(122, newReadCigarInfo.CorePositionEnd);
        assertEquals(32, newReadCigarInfo.FlankIndexEnd);
        assertEquals(132, newReadCigarInfo.FlankPositionEnd);
        assertEquals("29M", cigarElementsToStr(newReadCigarInfo.Cigar));

        // test 3: no HP and core end bases match
        // pos / index         10             20             30        40
        //           0123456789012345     67890     12345678901234567890123456789
        refBases =  "AACCGGTTAACCGATC" + "ACTCA" + "CTAGACGGTTAACCGGTT";
        readBases = "AACCGGTTAACCGATT" + "ACGCA" + "TTAGACGGTTAACCGGTT";
        //                 FFFFFFFFFF               FFFFFFFFFF

        refSequence = new RefSequence(100, refBases.getBytes());

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertTrue(checkEqual(readCigarInfo, newReadCigarInfo));
    }

    @Test
    public void testCoreExtensionSoftClips()
    {
        // test 1: left flank is in soft-clip and requires extension, making it invalid

        // pos / index                10             20             30
        //                  0123456789012345     67890     123456789012345678
        String refBases =  "AACCGGTTAACCGGTT" + "TACAT" + "TTAACCGGTTAACCGGTT";
        String readBases = refBases.substring(0, 18) + "G" + refBases.substring(19);

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

        // test 2: and now on the right in the core
        readAlignmentStart = 100;

        readCigarInfo = new ReadCigarInfo(
                readAlignmentStart, cigarElements, 106, 130, 116, 120,
                6, 30);

        cigarElements = Lists.newArrayList(
                new CigarElement(22, M),
                new CigarElement(17, S));

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNull(newReadCigarInfo);

        // test 3: accept truncated flank in append mode if we extend from M into S
        cigarElements = Lists.newArrayList(
                new CigarElement(31, M),
                new CigarElement(8, S));

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, true);

        assertNotNull(newReadCigarInfo);
        assertEquals(3, newReadCigarInfo.FlankIndexStart);
        assertEquals(103, newReadCigarInfo.FlankPositionStart);
        assertEquals(113, newReadCigarInfo.CorePositionStart);
        assertEquals(123, newReadCigarInfo.CorePositionEnd);
        assertEquals(30, newReadCigarInfo.FlankIndexEnd);
        assertEquals(130, newReadCigarInfo.FlankPositionEnd);
        assertEquals("28M", cigarElementsToStr(newReadCigarInfo.Cigar));

        // test 4: otherwise reject read context
        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);
        assertNull(newReadCigarInfo);


        // test 5: if core starts in soft-clip but does not require extension, still check/extend the other side
        refBases =  "AACCGGTTAACCGGTT" + "TACAT" + "ATAACCGGTTAACCGGTT";
        readBases = refBases.substring(0, 18) + "G" + refBases.substring(19);

        refSequence = new RefSequence(100, refBases.getBytes());

        cigarElements = Lists.newArrayList(
                new CigarElement(29, M),
                new CigarElement(10, S));

        newReadCigarInfo = extendUltimaCore(
                readBases.getBytes(), refSequence, readAlignmentStart, cigarElements, readCigarInfo, DEFAULT_FLANK_LENGTH, false);

        assertNotNull(newReadCigarInfo);
        assertEquals(3, newReadCigarInfo.FlankIndexStart);
        assertEquals(103, newReadCigarInfo.FlankPositionStart);
        assertEquals(113, newReadCigarInfo.CorePositionStart);
        assertEquals(120, newReadCigarInfo.CorePositionEnd);
        assertEquals(30, newReadCigarInfo.FlankIndexEnd);
        assertEquals(130, newReadCigarInfo.FlankPositionEnd);
        assertEquals("26M2S", cigarElementsToStr(newReadCigarInfo.Cigar));
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
        //                        0123456789012345678
        String middleRefBases =  "CTTTTTTTTTTAATTGCAA";
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
        int expectedCoreLength = 8;

        assertNotNull(readContext);
        assertEquals(10, readContext.CoreIndexStart);
        assertEquals(26, readContext.CoreIndexEnd);
        assertEquals("CTTTTTTTTTAATTGGC", readContext.coreStr());
        assertEquals("CTTTTTTTTTTAATTGC", readContext.refBases());
    }

    private static boolean checkEqual(final ReadCigarInfo first, final ReadCigarInfo second)
    {
        return first.FlankIndexStart == second.FlankIndexStart && first.FlankIndexEnd == second.FlankIndexEnd
                && first.CorePositionStart == second.CorePositionStart && first.CorePositionEnd == second.CorePositionEnd
                && first.FlankPositionStart == second.FlankPositionStart && first.FlankPositionEnd == second.FlankPositionEnd;
    }

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

}
