package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.redux.TestUtils.createFragmentCoords;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.redux.duplicate.FragmentCoords;

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

        FragmentCoords fragmentCoords = createFragmentCoords(read);

        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1:100_1:299:R_L", fragmentCoords.keyNonOriented());

        flipFirstInPair(read);

        fragmentCoords = createFragmentCoords(read);
        assertFalse(fragmentCoords.forwardFragment());

        // a mate upper read
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 100, true, false, null, false, TEST_READ_CIGAR);

        flipFirstInPair(read);

        fragmentCoords = createFragmentCoords(read);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertFalse(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1:100_1:299:R_U", fragmentCoords.keyNonOriented());

        // reverse orientation fragment with soft-clips at both ends
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "80M20S",
                CHR_1, 400, true, false, null, false, "10S80M");

        fragmentCoords = createFragmentCoords(read);
        assertEquals(199, fragmentCoords.PositionLower);
        assertEquals(390, fragmentCoords.PositionUpper);
        assertEquals(REVERSE, fragmentCoords.OrientLower);
        assertEquals(FORWARD, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1:199:R_1:390_L", fragmentCoords.keyNonOriented());

        // using supplementary data to get primary coords
        SupplementaryReadData suppAlignment = new SupplementaryReadData(
                CHR_1, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        read = createSamRecord(
                TEST_READ_ID, CHR_2, 1000, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, true, suppAlignment, true, TEST_READ_CIGAR);

        fragmentCoords = createFragmentCoords(read);
        assertNotNull(fragmentCoords.SuppReadInfo);
        assertEquals(1000, fragmentCoords.SuppReadInfo.UnclippedPosition);
        assertEquals(FORWARD, fragmentCoords.SuppReadInfo.Orient);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1:100_1:299:R_L_S", fragmentCoords.keyNonOriented());

        // test again with an unmapped mate
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                NO_CHROMOSOME_NAME, 100, false, false, null);

        fragmentCoords = createFragmentCoords(read);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(NO_POSITION, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1:100", fragmentCoords.keyNonOriented());

        read.setReadNegativeStrandFlag(true);
        fragmentCoords = createFragmentCoords(read);
        assertEquals(REVERSE, fragmentCoords.OrientLower);
        assertEquals("1:199:R", fragmentCoords.keyNonOriented());

        // the unmapped read
        read = createSamRecord(
                TEST_READ_ID, NO_CHROMOSOME_NAME, 100, TEST_READ_BASES, NO_CIGAR,
                CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
        read.setReadUnmappedFlag(true);

        fragmentCoords = createFragmentCoords(read);
        assertTrue(fragmentCoords.UnmappedSourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(NO_POSITION, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertFalse(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1:100_N", fragmentCoords.keyNonOriented());

        // an unpaired read
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                NO_CHROMOSOME_NAME, NO_POSITION, false, false, null, false, NO_CIGAR);
        read.setReadPairedFlag(false);

        fragmentCoords = createFragmentCoords(read);
        assertTrue(fragmentCoords.Unpaired);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(199, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.forwardFragment());
        assertEquals("1:100_199", fragmentCoords.keyNonOriented());

        // matching 5' positions, determine by first-in-pair
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);

        fragmentCoords = createFragmentCoords(read);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(100, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(FORWARD, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertEquals("1:100_1:100_L", fragmentCoords.keyNonOriented());

        flipFirstInPair(read);

        fragmentCoords = createFragmentCoords(read);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(100, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(FORWARD, fragmentCoords.OrientUpper);
        assertFalse(fragmentCoords.ReadIsLower);
        assertEquals("1:100_1:100_U", fragmentCoords.keyNonOriented());

    }

    @Test
    public void testFragmentCoordsUnpairedForward()
    {
        int readPos = 100;
        SAMRecord read = createSamRecordUnpaired(TEST_READ_ID, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, false, false, null);
        FragmentCoords fragmentCoords = createFragmentCoords(read);

        int positionLower = readPos;
        int positionUpper = readPos + TEST_READ_BASES.length() - 1;
        String expectedKey = format("%s:%d_%d", CHR_1, positionLower, positionUpper);

        assertEquals(CHR_1, fragmentCoords.ChromsomeLower);
        assertEquals(CHR_1, fragmentCoords.ChromsomeUpper);
        assertEquals(positionLower, fragmentCoords.PositionLower);
        assertEquals(positionUpper, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertEquals(FORWARD, fragmentCoords.FragmentOrient);
        assertTrue(fragmentCoords.ReadIsLower);
        assertNull(fragmentCoords.SuppReadInfo);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertTrue(fragmentCoords.Unpaired);
        assertEquals(expectedKey, fragmentCoords.Key);

        assertEquals(positionLower, fragmentCoords.readPosition());
        assertEquals(FORWARD, fragmentCoords.readOrientation());
    }

    @Test
    public void testFragmentCoordsUnpairedReverse()
    {
        int readPos = 100;
        SAMRecord read = createSamRecordUnpaired(TEST_READ_ID, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, true, false, null);
        FragmentCoords fragmentCoords = createFragmentCoords(read);

        int positionLower = readPos;
        int positionUpper = readPos + TEST_READ_BASES.length() - 1;
        String expectedKey = format("%s:%d_%d_R", CHR_1, positionLower, positionUpper);

        assertEquals(CHR_1, fragmentCoords.ChromsomeLower);
        assertEquals(CHR_1, fragmentCoords.ChromsomeUpper);
        assertEquals(positionLower, fragmentCoords.PositionLower);
        assertEquals(positionUpper, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertEquals(FORWARD, fragmentCoords.FragmentOrient);
        assertFalse(fragmentCoords.ReadIsLower);
        assertNull(fragmentCoords.SuppReadInfo);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertTrue(fragmentCoords.Unpaired);
        assertEquals(expectedKey, fragmentCoords.Key);

        assertEquals(positionUpper, fragmentCoords.readPosition());
        assertEquals(REVERSE, fragmentCoords.readOrientation());
    }

    @Test
    public void testFragmentCoordsUnpairedSuppForwardPrimaryForward()
    {
        String primaryChromosome = CHR_1;
        String suppChromosome = CHR_2;
        int primaryReadPos = 1_000;
        int suppReadPos = 2_000;

        SupplementaryReadData primarySuppReadData =
                new SupplementaryReadData(primaryChromosome, primaryReadPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 0);
        SAMRecord suppRead =
                createSamRecordUnpaired(TEST_READ_ID, suppChromosome, suppReadPos, TEST_READ_BASES, TEST_READ_CIGAR, false, true, primarySuppReadData);

        FragmentCoords fragmentCoords = createFragmentCoords(suppRead);

        int positionLower = primaryReadPos;
        int positionUpper = primaryReadPos + TEST_READ_BASES.length() - 1;
        String expectedKey = format("%s:%d_%d_S", primaryChromosome, positionLower, positionUpper);

        assertEquals(primaryChromosome, fragmentCoords.ChromsomeLower);
        assertEquals(primaryChromosome, fragmentCoords.ChromsomeUpper);
        assertEquals(positionLower, fragmentCoords.PositionLower);
        assertEquals(positionUpper, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertEquals(FORWARD, fragmentCoords.FragmentOrient);
        assertTrue(fragmentCoords.ReadIsLower);
        assertEquals(suppReadPos, fragmentCoords.SuppReadInfo.UnclippedPosition);
        assertEquals(FORWARD, fragmentCoords.SuppReadInfo.Orient);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertTrue(fragmentCoords.Unpaired);
        assertEquals(expectedKey, fragmentCoords.Key);

        assertEquals(suppReadPos, fragmentCoords.readPosition());
        assertEquals(FORWARD, fragmentCoords.readOrientation());
    }

    @Test
    public void testFragmentCoordsUnpairedSuppForwardPrimaryReverse()
    {
        String primaryChromosome = CHR_1;
        String suppChromosome = CHR_2;
        int primaryReadPos = 1_000;
        int suppReadPos = 2_000;

        SupplementaryReadData primarySuppReadData =
                new SupplementaryReadData(primaryChromosome, primaryReadPos, SUPP_NEG_STRAND, TEST_READ_CIGAR, 0);
        SAMRecord suppRead =
                createSamRecordUnpaired(TEST_READ_ID, suppChromosome, suppReadPos, TEST_READ_BASES, TEST_READ_CIGAR, false, true, primarySuppReadData);

        FragmentCoords fragmentCoords = createFragmentCoords(suppRead);

        int positionLower = primaryReadPos;
        int positionUpper = primaryReadPos + TEST_READ_BASES.length() - 1;
        String expectedKey = format("%s:%d_%d_R_S", primaryChromosome, positionLower, positionUpper);

        assertEquals(primaryChromosome, fragmentCoords.ChromsomeLower);
        assertEquals(primaryChromosome, fragmentCoords.ChromsomeUpper);
        assertEquals(positionLower, fragmentCoords.PositionLower);
        assertEquals(positionUpper, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertEquals(FORWARD, fragmentCoords.FragmentOrient);
        assertFalse(fragmentCoords.ReadIsLower);
        assertEquals(suppReadPos, fragmentCoords.SuppReadInfo.UnclippedPosition);
        assertEquals(FORWARD, fragmentCoords.SuppReadInfo.Orient);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertTrue(fragmentCoords.Unpaired);
        assertEquals(expectedKey, fragmentCoords.Key);

        assertEquals(suppReadPos, fragmentCoords.readPosition());
        assertEquals(FORWARD, fragmentCoords.readOrientation());
    }

    @Test
    public void testFragmentCoordsUnpairedSuppReversePrimaryForward()
    {
        String primaryChromosome = CHR_1;
        String suppChromosome = CHR_2;
        int primaryReadPos = 1_000;
        int suppReadPos = 2_000;

        SupplementaryReadData primarySuppReadData =
                new SupplementaryReadData(primaryChromosome, primaryReadPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 0);
        SAMRecord suppRead =
                createSamRecordUnpaired(TEST_READ_ID, suppChromosome, suppReadPos, TEST_READ_BASES, TEST_READ_CIGAR, true, true, primarySuppReadData);

        FragmentCoords fragmentCoords = createFragmentCoords(suppRead);

        int positionLower = primaryReadPos;
        int positionUpper = primaryReadPos + TEST_READ_BASES.length() - 1;
        String expectedKey = format("%s:%d_%d_S", primaryChromosome, positionLower, positionUpper);

        assertEquals(primaryChromosome, fragmentCoords.ChromsomeLower);
        assertEquals(primaryChromosome, fragmentCoords.ChromsomeUpper);
        assertEquals(positionLower, fragmentCoords.PositionLower);
        assertEquals(positionUpper, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertEquals(FORWARD, fragmentCoords.FragmentOrient);
        assertTrue(fragmentCoords.ReadIsLower);
        assertEquals(suppReadPos + TEST_READ_BASES.length() - 1, fragmentCoords.SuppReadInfo.UnclippedPosition);
        assertEquals(REVERSE, fragmentCoords.SuppReadInfo.Orient);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertTrue(fragmentCoords.Unpaired);
        assertEquals(expectedKey, fragmentCoords.Key);

        assertEquals(suppReadPos + TEST_READ_BASES.length() - 1, fragmentCoords.readPosition());
        assertEquals(REVERSE, fragmentCoords.readOrientation());
    }

    @Test
    public void testFragmentCoordsUnpairedSuppReversePrimaryReverse()
    {
        String primaryChromosome = CHR_1;
        String suppChromosome = CHR_2;
        int primaryReadPos = 1_000;
        int suppReadPos = 2_000;

        SupplementaryReadData primarySuppReadData =
                new SupplementaryReadData(primaryChromosome, primaryReadPos, SUPP_NEG_STRAND, TEST_READ_CIGAR, 0);
        SAMRecord suppRead =
                createSamRecordUnpaired(TEST_READ_ID, suppChromosome, suppReadPos, TEST_READ_BASES, TEST_READ_CIGAR, true, true, primarySuppReadData);

        FragmentCoords fragmentCoords = createFragmentCoords(suppRead);

        int positionLower = primaryReadPos;
        int positionUpper = primaryReadPos + TEST_READ_BASES.length() - 1;
        String expectedKey = format("%s:%d_%d_R_S", primaryChromosome, positionLower, positionUpper);

        assertEquals(primaryChromosome, fragmentCoords.ChromsomeLower);
        assertEquals(primaryChromosome, fragmentCoords.ChromsomeUpper);
        assertEquals(positionLower, fragmentCoords.PositionLower);
        assertEquals(positionUpper, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertEquals(FORWARD, fragmentCoords.FragmentOrient);
        assertFalse(fragmentCoords.ReadIsLower);
        assertEquals(suppReadPos + TEST_READ_BASES.length() - 1, fragmentCoords.SuppReadInfo.UnclippedPosition);
        assertEquals(REVERSE, fragmentCoords.SuppReadInfo.Orient);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertTrue(fragmentCoords.Unpaired);
        assertEquals(expectedKey, fragmentCoords.Key);

        assertEquals(suppReadPos + TEST_READ_BASES.length() - 1, fragmentCoords.readPosition());
        assertEquals(REVERSE, fragmentCoords.readOrientation());
    }
}
