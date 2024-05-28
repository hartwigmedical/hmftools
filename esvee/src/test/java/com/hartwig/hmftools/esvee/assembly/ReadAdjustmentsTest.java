package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.AssemblyConstants.POLY_G_TRIM_LENGTH;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_CIGAR_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;

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
        return ReadAdjustments.convertEdgeIndelsToSoftClip(read, 6, 15);
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
        assertEquals(100, read.unclippedStart());
        assertEquals(134, read.unclippedEnd());
        assertEquals(118, read.indelImpliedAlignmentStart());
        assertEquals(117, read.indelImpliedAlignmentEnd());
        assertEquals(94, read.minUnclippedStart());
        assertEquals(140, read.maxUnclippedEnd());

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
        assertEquals(100, read.unclippedStart());
        assertEquals(130, read.unclippedEnd());
        assertEquals(116, read.indelImpliedAlignmentStart());
        assertEquals(109, read.indelImpliedAlignmentEnd());
        assertEquals(100, read.minUnclippedStart());
        assertEquals(130, read.maxUnclippedEnd());
    }
}