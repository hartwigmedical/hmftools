package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_CIGAR_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.hasLineTail;
import static com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments.trimAdapterBases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.prep.ReadFilters;

import org.junit.Test;

public class ReadAdjustmentsTest
{
    @Test
    public void testPolyGTrimming()
    {
        String polyGSection = "GGGGGGGG";
        String polyCSection = "CCCCCCCC";
        String otherBases = "AAAACCCCGAGATTTT";

        // wrong end
        String readBases = polyGSection + otherBases;
        Read read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 0));
        assertFalse(ReadAdjustments.trimPolyGSequences(read));

        // too few Gs
        readBases = otherBases + polyGSection.substring(0, 3);
        read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 0));
        assertFalse(ReadAdjustments.trimPolyGSequences(read));

        // trimmed on the 3' end
        readBases = otherBases + polyGSection.substring(0, POLY_G_TRIM_LENGTH);
        read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 0));
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(otherBases.length(), read.basesLength());
        assertEquals(read.alignmentStart() + otherBases.length() - 1, read.alignmentEnd());
        assertEquals(otherBases, read.getBasesString());
        assertEquals(makeCigarString(otherBases, 0, 0), read.cigarString());

        // again with a soft-clipping at the end
        readBases = otherBases + polyGSection.substring(0, 8);
        read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 4));
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(otherBases.length(), read.basesLength());
        assertEquals(read.alignmentStart() + otherBases.length() - 1, read.alignmentEnd());
        assertEquals(otherBases, read.getBasesString());
        assertEquals(makeCigarString(otherBases, 0, 0), read.cigarString());

        // negative strand - trimmed at the start
        readBases = polyCSection.substring(0, 6) + otherBases;
        read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 3, 0));
        read.bamRecord().setReadNegativeStrandFlag(true);
        assertEquals(97, read.unclippedStart());

        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(103, read.alignmentStart());
        assertEquals(103, read.unclippedStart());
        assertEquals(103 + otherBases.length() - 1, read.alignmentEnd());
        assertEquals(otherBases, read.getBasesString());
        assertEquals(makeCigarString(otherBases, 0, 0), read.cigarString());

        // test trimming into delete
        readBases = otherBases + polyGSection.substring(0, 5);
        read = createRead(TEST_READ_ID, 100, readBases, "16M2D5M");
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(115, read.alignmentEnd());
        assertEquals("16M", read.cigarString());

        readBases = otherBases + polyGSection.substring(0, 6);
        read = createRead(TEST_READ_ID, 100, readBases, "16M2D5M");
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(114, read.alignmentEnd());
        assertEquals("15M", read.cigarString());

        // again from the other end
        readBases = polyCSection.substring(0, 5) + otherBases;
        read = createRead(TEST_READ_ID, 100, readBases, "5M2D16M");
        read.bamRecord().setReadNegativeStrandFlag(true);
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(107, read.alignmentStart());
        assertEquals("16M", read.cigarString());

        readBases = polyCSection.substring(0, 6) + otherBases;
        read = createRead(TEST_READ_ID, 100, readBases, "5M2D16M");
        read.bamRecord().setReadNegativeStrandFlag(true);
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(108, read.alignmentStart());
        assertEquals("15M", read.cigarString());
    }

    @Test
    public void testLowBaseFiltering()
    {
        String readBases = REF_BASES_RANDOM_100;
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());

        byte lowQualBase = 10;
        for(int i = 0; i < 50; ++i)
        {
            baseQualities[i] = lowQualBase;
        }

        Read read = createRead(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);

        assertFalse(ReadFilters.filterLowQualRead(read.bamRecord()));

        baseQualities[50] = lowQualBase;
        assertTrue(ReadFilters.filterLowQualRead(read.bamRecord()));
    }

    @Test
    public void testSoftClipLowBaseQualTrimming()
    {
        // wrong end
        String readBases = REF_BASES_RANDOM_100;

        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());

        int trimLength = 8;
        byte lowQualBase = 10;
        for(int i = 0; i < trimLength; ++i)
        {
            baseQualities[i] = lowQualBase;
            baseQualities[baseQualities.length - i - 1] = lowQualBase;
        }

        Read read = createRead(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);

        assertFalse(ReadAdjustments.trimLowQualSoftClipBases(read)); // nothing without soft-clips

        read = createRead(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualSoftClipBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(189, read.alignmentEnd());
        assertEquals(191, read.unclippedEnd());
        assertEquals("10S80M2S", read.cigarString());

        read = createRead(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));

        for(int i = 90; i < baseQualities.length; ++i)
        {
            baseQualities[i] = DEFAULT_BASE_QUAL;
        }

        baseQualities[95] = lowQualBase;
        baseQualities[98] = lowQualBase;
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualSoftClipBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(189, read.alignmentEnd());
        assertEquals(197, read.unclippedEnd());
        assertEquals("10S80M8S", read.cigarString());

        read = createRead(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));
        read.bamRecord().setReadNegativeStrandFlag(true);
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualSoftClipBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(108, read.unclippedStart());
        assertEquals("2S80M10S", read.cigarString());
    }

    @Test
    public void testLowBaseQualTrimming()
    {
        String readBases = REF_BASES_RANDOM_100;

        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());

        // low qual will be 70% of outer bases
        byte lowQualBase = 10;
        for(int i = 3; i < 10; ++i)
        {
            baseQualities[i] = lowQualBase;
            baseQualities[baseQualities.length - i - 1] = lowQualBase;
        }

        Read read = createRead(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);

        read.trimLowQualBases();

        assertTrue(read.lowQualTrimmed());
        assertEquals("90M", read.cigarString());
        assertEquals(10, read.baseTrimCount());
        assertEquals(189, read.alignmentEnd());

        read = createRead(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);
        read.bamRecord().setReadNegativeStrandFlag(true);

        read.trimLowQualBases();

        assertTrue(read.lowQualTrimmed());
        assertEquals("90M", read.cigarString());
        assertEquals(10, read.baseTrimCount());
        assertEquals(110, read.alignmentStart());
    }

    @Test
    public void testAdapter()
    {
        String readBases = REF_BASES_RANDOM_100;

        Read read = createRead(TEST_READ_ID, 110, readBases, "20S60M20S");
        read.bamRecord().setMateNegativeStrandFlag(true);

        Read mateRead = createRead(TEST_READ_ID, 100, readBases, "20S60M20S");
        mateRead.bamRecord().setReadNegativeStrandFlag(true);

        read.setMateRead(mateRead);
        mateRead.setMateRead(read);

        trimAdapterBases(read);
        trimAdapterBases(mateRead);

        assertEquals(10, read.baseTrimCount());
        assertEquals("20S60M10S", read.cigarString());
        assertEquals(179, read.unclippedEnd());

        assertEquals(10, mateRead.baseTrimCount());
        assertEquals("10S60M20S", mateRead.cigarString());
        assertEquals(90, mateRead.unclippedStart());
    }

    @Test
    public void testMarkLineTailRead()
    {
        String lineSequence = "AAAAAAAAAAAAAAAA";
        String nonLineSequence = "GCTGCTGTCGTGTCC";
        String softClipBases = nonLineSequence + lineSequence;
        String readBases = softClipBases + REF_BASES_RANDOM_100;

        assertTrue(hasLineTail(readBases.getBytes(), softClipBases.length() - 1, true, LINE_BASE_A));

        lineSequence = "AAAAAAGAAAAAATAAAA";
        softClipBases = nonLineSequence + lineSequence;
        readBases = softClipBases + REF_BASES_RANDOM_100;

        assertTrue(hasLineTail(readBases.getBytes(), softClipBases.length() - 1, true, LINE_BASE_A));

        // too many non-line bases
        lineSequence = "AAAAAAGAAACAATAAAA";
        softClipBases = nonLineSequence + lineSequence;
        readBases = softClipBases + REF_BASES_RANDOM_100;

        assertFalse(hasLineTail(readBases.getBytes(), softClipBases.length() - 1, true, LINE_BASE_A));

        // an extension of a local repeat
        lineSequence = "AAAAAAAAAAAAAAAA";
        softClipBases = nonLineSequence + lineSequence;
        readBases = softClipBases + "AAAAAAAAAA" + REF_BASES_RANDOM_100;

        assertFalse(hasLineTail(readBases.getBytes(), softClipBases.length() - 1, true, LINE_BASE_A));

        // reverse side
        lineSequence = "TTTTTCTTTTTTTGTTTT";
        softClipBases = lineSequence + nonLineSequence;
        readBases = REF_BASES_RANDOM_100 + softClipBases;

        assertTrue(hasLineTail(readBases.getBytes(), REF_BASES_RANDOM_100.length(), false, LINE_BASE_T));

        readBases = REF_BASES_RANDOM_100 + "TTTTTTTTTT" + softClipBases;

        assertFalse(hasLineTail(readBases.getBytes(), REF_BASES_RANDOM_100.length() + 10, false, LINE_BASE_T));
    }

    @Test
    public void testLineLowQualTrimming()
    {
        String lineSequence = "AAAAAAAAAAAAAAAA";
        String softClipBases = "GCTGCTGTCGTGTCC" + lineSequence;
        String readBases = softClipBases + REF_BASES_RANDOM_100;

        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());

        int lowQualLength = readBases.length() - 2;
        byte lowQualBase = 10;
        for(int i = 0; i < lowQualLength; ++i)
        {
            baseQualities[i] = lowQualBase;
        }

        Read read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, softClipBases.length(), 0));
        read.bamRecord().setBaseQualities(baseQualities);
        read.bamRecord().setReadNegativeStrandFlag(true);

        ReadAdjustments.markLineSoftClips(read);
        assertTrue(ReadAdjustments.trimLowQualSoftClipBases(read));
        assertEquals(84, read.unclippedStart());
        assertEquals("16S100M", read.cigarString());

        // keep 2 extra bases since the non-line bases are further in
        softClipBases = "GCTGCTGTCGTGT" + lineSequence + "CC";
        readBases = softClipBases + REF_BASES_RANDOM_100;

        read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, softClipBases.length(), 0));
        read.bamRecord().setBaseQualities(baseQualities);
        read.bamRecord().setReadNegativeStrandFlag(true);

        ReadAdjustments.markLineSoftClips(read);
        assertTrue(ReadAdjustments.trimLowQualSoftClipBases(read));
        assertEquals(82, read.unclippedStart());
        assertEquals("18S100M", read.cigarString());

        // test no LINE trimming if an extension of a matching repeat
        readBases = softClipBases + lineSequence + REF_BASES_RANDOM_100;

        read = createRead(TEST_READ_ID, 100, readBases, makeCigarString(readBases, softClipBases.length(), 0));
        read.bamRecord().setBaseQualities(baseQualities);
        read.bamRecord().setReadNegativeStrandFlag(true);

        ReadAdjustments.markLineSoftClips(read);
        assertFalse(read.hasLineTail());


        // test the other orientation
        lineSequence = Nucleotides.reverseComplementBases(lineSequence);
        softClipBases = "CC" + lineSequence + "GCTGCTGTCGTGTCC";
        readBases = REF_BASES_RANDOM_100 + softClipBases;

        baseQualities = buildDefaultBaseQuals(readBases.length());

        for(int i = 0; i < readBases.length(); ++i)
        {
            baseQualities[i] = lowQualBase;
        }

        read = createRead(TEST_READ_ID, 101, readBases, makeCigarString(readBases, 0, softClipBases.length()));
        read.bamRecord().setBaseQualities(baseQualities);

        ReadAdjustments.markLineSoftClips(read);
        assertTrue(ReadAdjustments.trimLowQualSoftClipBases(read));
        assertEquals(200, read.alignmentEnd());
        assertEquals(218, read.unclippedEnd());
        assertEquals("100M18S", read.cigarString());


        // expansion of local repeat
        readBases = REF_BASES_RANDOM_100 + lineSequence + softClipBases;
        read = createRead(TEST_READ_ID, 101, readBases, makeCigarString(readBases, 0, softClipBases.length()));
        read.bamRecord().setBaseQualities(baseQualities);

        ReadAdjustments.markLineSoftClips(read);
        assertFalse(read.hasLineTail());
    }
}