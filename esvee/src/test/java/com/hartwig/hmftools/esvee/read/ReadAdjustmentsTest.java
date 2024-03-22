package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_CIGAR_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.esvee.AssemblyConstants;

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
        readBases = otherBases + polyGSection.substring(0, AssemblyConstants.POLY_G_TRIM_LENGTH);
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
    }

    @Test
    public void testLowBaseQualTrimming()
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

        assertFalse(ReadAdjustments.trimLowQualBases(read)); // nothing without soft-clips

        read = createRead(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(189, read.alignmentEnd());
        assertEquals(191, read.unclippedEnd());
        assertEquals("10S80M2S", read.cigarString());

        // hits > 30% at index 5 but not again
        read = createRead(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));

        for(int i = 90; i < baseQualities.length; ++i)
        {
            baseQualities[i] = DEFAULT_BASE_QUAL;
        }

        baseQualities[95] = lowQualBase;
        baseQualities[98] = lowQualBase;
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(189, read.alignmentEnd());
        assertEquals(194, read.unclippedEnd());
        assertEquals("10S80M5S", read.cigarString());

        read = createRead(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));
        read.bamRecord().setReadNegativeStrandFlag(true);
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(108, read.unclippedStart());
        assertEquals("2S80M10S", read.cigarString());
    }

    private static boolean convertEdgeIndelsToSoftClip(final Read read, boolean allowDoubleConversion)
    {
        return ReadAdjustments.convertEdgeIndelsToSoftClip(read, 6, 15, allowDoubleConversion);
    }

    @Test
    public void testImpliedEdgeIndelsToSoftClip()
    {
        // the cigar remains unch but the implied new alignments are calculated and the unclipped alignments stored

        // converts both sides of the insert
        String cigar = "18M6I17M";
        String readBases = REF_BASES_RANDOM_100.substring(0, 40);
        Read read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertTrue(convertEdgeIndelsToSoftClip(read, true));
        assertEquals(100, read.alignmentStart());
        assertEquals(134, read.alignmentEnd());
        assertEquals(94, read.unclippedStart());
        assertEquals(140, read.unclippedEnd());
        assertEquals(118, read.indelImpliedAlignmentStart());
        assertEquals(117, read.indelImpliedAlignmentEnd());

        // indel too short
        cigar = "10M5I10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 25);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertFalse(convertEdgeIndelsToSoftClip(read, true));

        // converts both sides of the delete
        cigar = "10M6D15M";
        readBases = REF_BASES_RANDOM_100.substring(0, 31);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(130, read.alignmentEnd());
        assertTrue(convertEdgeIndelsToSoftClip(read, true));
        assertEquals(100, read.alignmentStart());
        assertEquals(130, read.alignmentEnd());
        assertEquals(106, read.unclippedStart());
        assertEquals(124, read.unclippedEnd());
        assertEquals(116, read.indelImpliedAlignmentStart());
        assertEquals(109, read.indelImpliedAlignmentEnd());
    }

    @Test
    public void testEdgeIndelsToSoftClip()
    {
        String cigar = "18M6I17M";
        String readBases = REF_BASES_RANDOM_100.substring(0, 40);
        Read read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertTrue(convertEdgeIndelsToSoftClip(read, false));
        assertEquals(117, read.alignmentEnd());
        assertEquals(140, read.unclippedEnd());

        // indel too short
        cigar = "10M5I10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 25);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertFalse(convertEdgeIndelsToSoftClip(read, false));

        // convert a left-edge delete
        cigar = "10M6D15M";
        readBases = REF_BASES_RANDOM_100.substring(0, 31);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(100, read.unclippedStart());
        assertEquals(130, read.alignmentEnd());
        assertTrue(convertEdgeIndelsToSoftClip(read, false));
        assertEquals(116, read.alignmentStart());
        assertEquals(106, read.unclippedStart());
        assertEquals(130, read.alignmentEnd());
        assertEquals("10S15M", read.cigarString());

        // convert a right-edge delete
        cigar = "15M6D10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 31);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(130, read.alignmentEnd());
        assertTrue(convertEdgeIndelsToSoftClip(read, false));
        assertEquals(100, read.alignmentStart());
        assertEquals(114, read.alignmentEnd());
        assertEquals(124, read.unclippedEnd());
        assertEquals("15M10S", read.cigarString());

        // convert a left-edge insert
        cigar = "10M6I15M";
        readBases = REF_BASES_RANDOM_100.substring(0, 31);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(124, read.alignmentEnd());
        assertEquals(100, read.unclippedStart());

        assertTrue(convertEdgeIndelsToSoftClip(read, false));

        assertEquals(110, read.alignmentStart());
        assertEquals(94, read.unclippedStart());
        assertEquals(124, read.alignmentEnd());
        assertEquals("16S15M", read.cigarString());

        // right edge
        cigar = "20M6I15M";
        readBases = REF_BASES_RANDOM_100.substring(0, 41);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(134, read.alignmentEnd());
        assertEquals(134, read.unclippedEnd());

        assertTrue(convertEdgeIndelsToSoftClip(read, false));

        assertEquals(119, read.alignmentEnd());
        assertEquals(140, read.unclippedEnd());
        assertEquals("20M21S", read.cigarString());

        // only the shorter side will be converted
        cigar = "12M8D10M6I10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 36);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(139, read.alignmentEnd());
        assertEquals(139, read.unclippedEnd());

        assertTrue(convertEdgeIndelsToSoftClip(read, false));

        assertEquals(100, read.alignmentStart());
        assertEquals(100, read.unclippedStart());
        assertEquals(129, read.alignmentEnd());
        assertEquals(145, read.unclippedEnd());
        assertEquals("12M8D10M16S", read.cigarString());
    }
}