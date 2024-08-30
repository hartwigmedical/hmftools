// TODO: REVIEW
package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.sage.SageConfig;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaCoreExtenderTest
{
    @Test
    public void testBuildContextNoCoreExtension()
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
        String cigar = String.valueOf(readBases.length()) + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readBases, cigar, false, false, null);

        int varIndexInRead = readPaddingSize + variantMiddleIndex;

        RefSequence refSequence = new RefSequence(1, refBases.getBytes());

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA);
        VariantReadContextBuilder readContextBuilder = new VariantReadContextBuilder(config, DEFAULT_FLANK_LENGTH);

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
    public void testBuildContextWithExtension()
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
        String cigar = String.valueOf(readBases.length()) + "M";
        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readBases, cigar, false, false, null);

        int varIndexInRead = readPaddingSize + variantMiddleIndex;

        RefSequence refSequence = new RefSequence(1, refBases.getBytes());

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA);
        VariantReadContextBuilder readContextBuilder = new VariantReadContextBuilder(config, DEFAULT_FLANK_LENGTH);

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
