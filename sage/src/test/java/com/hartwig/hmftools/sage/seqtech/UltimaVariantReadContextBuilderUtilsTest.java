package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.sage.common.TestUtils.setUltimaSequencing;
import static com.hartwig.hmftools.sage.common.UltimaVariantReadContextBuilderUtils.isAdjacentToLongHomopolymer;
import static com.hartwig.hmftools.sage.common.UltimaVariantReadContextBuilderUtils.isMsiIndelOfType;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaVariantReadContextBuilderUtilsTest
{
    public UltimaVariantReadContextBuilderUtilsTest()
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
    public void testBuildContextUltimaNoCoreExtension()
    {
        String middleRefBases = "CTCGGCCTCCCGGAGGTGCCGG";
        String middleReadBases = "CTCGGCCTCCCTGAGGTGCCGG";
        int variantMiddleIndex = 11;
        int refPaddingSize = 1000;
        int readPaddingSize = 100;
        String refPadding = "A".repeat(refPaddingSize);
        String readPadding = "A".repeat(readPaddingSize);
        String refBases = refPadding + middleRefBases + refPadding;
        String readBases = readPadding + middleReadBases + readPadding;

        SimpleVariant variant = new SimpleVariant(CHR_1, 1 + refPaddingSize + variantMiddleIndex, "G", "T");

        int alignmentStart = 1 + refPaddingSize - readPaddingSize;
        String cigar = readBases.length() + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readBases, cigar, false, false, null);

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
        int refPaddingSize = 1000;
        int readPaddingSize = 100;
        String refBases = "A".repeat(refPaddingSize) + middleRefBases + "T".repeat(refPaddingSize);
        String readBases = "A".repeat(readPaddingSize) + middleReadBases + "T".repeat(readPaddingSize);

        SimpleVariant variant = new SimpleVariant(CHR_1, 1 + refPaddingSize + variantMiddleIndex, "T", "G");

        int alignmentStart = 1 + refPaddingSize - readPaddingSize;
        String cigar = readBases.length() + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readBases, cigar, false, false, null);

        int varIndexInRead = readPaddingSize + variantMiddleIndex;

        RefSequence refSequence = new RefSequence(1, refBases.getBytes());

        // SageConfig config = new SageConfig(false, SequencingType.ULTIMA);
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
