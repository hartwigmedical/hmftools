package com.hartwig.hmftools.sage.ultima;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.sage.common.TestUtils.setUltimaSequencing;
import static com.hartwig.hmftools.sage.common.UltimaVariantReadContextBuilderUtils.isAdjacentToLongHomopolymer;
import static com.hartwig.hmftools.sage.common.UltimaVariantReadContextBuilderUtils.isMsiIndelOfType;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import org.junit.After;
import org.junit.Ignore;
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
    public void testIsMsiIndelOfTypeNotIndel()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(false);

        assertFalse(isMsiIndelOfType(mockVariant, null));
    }

    @Test
    public void testIsMsiIndelOfTypeNoUnits()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);
        when(mockVariant.isInsert()).thenReturn(true);
        when(mockVariant.ref()).thenReturn("A");
        when(mockVariant.alt()).thenReturn("AA");

        List<String> units = Lists.newArrayList();
        assertFalse(isMsiIndelOfType(mockVariant, units));
    }

    @Test
    public void testIsMsiIndelOfTypeSingleUnitNotMsiIndel()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);
        when(mockVariant.isInsert()).thenReturn(true);
        when(mockVariant.ref()).thenReturn("A");
        when(mockVariant.alt()).thenReturn("AAB");

        List<String> units = Lists.newArrayList("A");
        assertFalse(isMsiIndelOfType(mockVariant, units));
    }

    @Test
    public void testIsMsiIndelOfTypeSingleUnitIsMsiIndel()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);
        when(mockVariant.isInsert()).thenReturn(true);
        when(mockVariant.ref()).thenReturn("A");
        when(mockVariant.alt()).thenReturn("AABABAB");

        List<String> units = Lists.newArrayList("AB");
        assertTrue(isMsiIndelOfType(mockVariant, units));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerNotLongEnough()
    {
        int longLength = 10;
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn("A".repeat(2*longLength - 1).getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength - 1, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerOnLeft()
    {
        int longLength = 10;
        String readBases = "A".repeat(longLength) + "C" + "T".repeat(longLength - 1);
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertTrue(isAdjacentToLongHomopolymer(mockRead, longLength, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerOnRight()
    {
        int longLength = 10;
        String readBases = "T".repeat(longLength - 1) + "C" + "A".repeat(longLength);
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertTrue(isAdjacentToLongHomopolymer(mockRead, longLength - 1, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerNotOnLeft()
    {
        int longLength = 10;
        String readBases = "C" + "A".repeat(2*longLength - 1);
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerNotOnRight()
    {
        int longLength = 10;
        String readBases = "A".repeat(2*longLength - 1) + "C";
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength - 1, longLength));
    }

    @Ignore
    @Test
    public void testIsAdjacentToLongHomopolymerNHomopolymerOnLeft()
    {
        int longLength = 10;
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn("N".repeat(2*longLength).getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength, longLength));
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
