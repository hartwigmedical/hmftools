package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_CIGAR_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.SvConstants;

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
        Read read = createSamRecord(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 0));
        assertFalse(ReadAdjustments.trimPolyGSequences(read));

        // too few Gs
        readBases = otherBases + polyGSection.substring(0, 3);
        read = createSamRecord(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 0));
        assertFalse(ReadAdjustments.trimPolyGSequences(read));

        // trimmed on the 3' end
        readBases = otherBases + polyGSection.substring(0, SvConstants.POLY_G_TRIM_LENGTH);
        read = createSamRecord(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 0));
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(otherBases.length(), read.basesLength());
        assertEquals(read.alignmentStart() + otherBases.length() - 1, read.alignmentEnd());
        assertEquals(otherBases, read.getBasesString());
        assertEquals(makeCigarString(otherBases, 0, 0), read.cigarString());

        // again with a soft-clipping at the end
        readBases = otherBases + polyGSection.substring(0, 8);
        read = createSamRecord(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 0, 4));
        assertTrue(ReadAdjustments.trimPolyGSequences(read));

        assertEquals(otherBases.length(), read.basesLength());
        assertEquals(read.alignmentStart() + otherBases.length() - 1, read.alignmentEnd());
        assertEquals(otherBases, read.getBasesString());
        assertEquals(makeCigarString(otherBases, 0, 0), read.cigarString());

        // negative strand - trimmed at the start
        readBases = polyCSection.substring(0, 6) + otherBases;
        read = createSamRecord(TEST_READ_ID, 100, readBases, makeCigarString(readBases, 3, 0));
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

        Read read = createSamRecord(TEST_READ_ID, 100, readBases, TEST_CIGAR_100);
        read.bamRecord().setBaseQualities(baseQualities);

        assertFalse(ReadAdjustments.trimLowQualBases(read)); // nothing without soft-clips

        read = createSamRecord(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(189, read.alignmentEnd());
        assertEquals(191, read.unclippedEnd());
        assertEquals("10S80M2S", read.cigarString());

        // hits > 30% at index 5 but not again
        read = createSamRecord(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));

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

        read = createSamRecord(TEST_READ_ID, 110, readBases, makeCigarString(readBases, 10, 10));
        read.bamRecord().setReadNegativeStrandFlag(true);
        read.bamRecord().setBaseQualities(baseQualities);

        assertTrue(ReadAdjustments.trimLowQualBases(read));
        assertEquals(110, read.alignmentStart());
        assertEquals(108, read.unclippedStart());
        assertEquals("2S80M10S", read.cigarString());
    }

    @Test
    public void testEdgeIndelsToSoftClip()
    {
        // to far from edge
        String cigar = "17M6I17M";
        String readBases = REF_BASES_RANDOM_100.substring(0, 40);
        Read read = createSamRecord(TEST_READ_ID, 100, readBases, cigar);
        assertFalse(ReadAdjustments.convertEdgeIndelsToSoftClip(read));

        // indel too short
        cigar = "10M5I10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 25);
        read = createSamRecord(TEST_READ_ID, 100, readBases, cigar);
        assertFalse(ReadAdjustments.convertEdgeIndelsToSoftClip(read));

        // convert a left-edge indel
        cigar = "10M6I15M";
        readBases = REF_BASES_RANDOM_100.substring(0, 31);
        read = createSamRecord(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(124, read.alignmentEnd());
        assertEquals(100, read.unclippedStart());

        assertTrue(ReadAdjustments.convertEdgeIndelsToSoftClip(read));

        assertEquals(110, read.alignmentStart());
        assertEquals(94, read.unclippedStart());
        assertEquals(124, read.alignmentEnd());
        assertEquals("16S15M", read.cigarString());

        // right edge
        cigar = "20M6I15M";
        readBases = REF_BASES_RANDOM_100.substring(0, 41);
        read = createSamRecord(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(134, read.alignmentEnd());
        assertEquals(134, read.unclippedEnd());

        assertTrue(ReadAdjustments.convertEdgeIndelsToSoftClip(read));

        assertEquals(119, read.alignmentEnd());
        assertEquals(140, read.unclippedEnd());
        assertEquals("20M21S", read.cigarString());

        // convert both at once
        cigar = "10M8D10M6I10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 36);
        read = createSamRecord(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(137, read.alignmentEnd());
        assertEquals(137, read.unclippedEnd());

        assertTrue(ReadAdjustments.convertEdgeIndelsToSoftClip(read));

        assertEquals(110, read.alignmentStart());
        assertEquals(100, read.unclippedStart());
        assertEquals(119, read.alignmentEnd());
        assertEquals(135, read.unclippedEnd());
        assertEquals("10S10M16S", read.cigarString());


    }
}