package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.redux.common.FragmentCoords;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class FragmentUtilsTest
{
    @Test
    public void testFragmentCoords()
    {
        SAMRecord read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, false, null, true, TEST_READ_CIGAR);

        FragmentCoords fragmentCoords = FragmentCoords.fromRead(read);

        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1_100_1_299_R", fragmentCoords.Key);

        flipFirstInPair(read);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertFalse(fragmentCoords.forwardFragment());

        // a mate upper read
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 100, true, false, null, false, TEST_READ_CIGAR);

        flipFirstInPair(read);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertFalse(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1_100_1_299_R", fragmentCoords.Key);

        // reverse orientation fragment with soft-clips at both ends
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "80M20S",
                CHR_1, 400, true, false, null, false, "10S80M");

        fragmentCoords = FragmentCoords.fromRead(read);
        assertEquals(199, fragmentCoords.PositionLower);
        assertEquals(390, fragmentCoords.PositionUpper);
        assertEquals(REVERSE, fragmentCoords.OrientLower);
        assertEquals(FORWARD, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1_199_R_1_390", fragmentCoords.Key);

        // using supplementary data to get primary coords
        SupplementaryReadData suppAlignment = new SupplementaryReadData(
                CHR_1, 100, SupplementaryReadData.SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        read = createSamRecord(
                TEST_READ_ID, CHR_2, 1000, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, true, suppAlignment, true, TEST_READ_CIGAR);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertTrue(fragmentCoords.SupplementarySourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1_100_1_299_R", fragmentCoords.Key);

        // test again with an unmapped mate
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                NO_CHROMOSOME_NAME, 100, false, false, null);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(NO_POSITION, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1_100", fragmentCoords.Key);

        read.setReadNegativeStrandFlag(true);
        fragmentCoords = FragmentCoords.fromRead(read);
        assertEquals(REVERSE, fragmentCoords.OrientLower);
        assertEquals("1_199_R", fragmentCoords.Key);

        // the unmapped read
        read = createSamRecord(
                TEST_READ_ID, NO_CHROMOSOME_NAME, 100, TEST_READ_BASES, NO_CIGAR,
                CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
        read.setReadUnmappedFlag(true);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertTrue(fragmentCoords.UnmappedSourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(NO_POSITION, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertFalse(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1_100_U", fragmentCoords.Key);
    }
}
