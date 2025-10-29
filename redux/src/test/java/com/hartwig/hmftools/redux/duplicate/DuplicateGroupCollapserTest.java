package com.hartwig.hmftools.redux.duplicate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.redux.duplicate.CollapseUtils.collapseToNonOrientedKeyWithoutCoordinates;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroupCollapserTest
{
    private static final String TEST_READ_BASES = "A".repeat(100);
    private static final String TEST_READ_CIGAR = "100M";

    @Test
    public void testCollapseToKeyWithoutCoordinatesUnmappedRead()
    {
        checkCollapseToKeyWithoutCoordinatesUnmappedRead(true, true);
        checkCollapseToKeyWithoutCoordinatesUnmappedRead(true, false);
        checkCollapseToKeyWithoutCoordinatesUnmappedRead(false, true);
        checkCollapseToKeyWithoutCoordinatesUnmappedRead(false, false);
    }

    private void checkCollapseToKeyWithoutCoordinatesUnmappedRead(boolean firstOfPair, boolean mateReversed)
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                "READ_001", NO_CHROMOSOME_NAME, NO_POSITION, TEST_READ_BASES, NO_CIGAR, CHR_1,
                1_000, false, false, null,
                mateReversed, TEST_READ_CIGAR);
        read.setReadUnmappedFlag(true);
        read.setFirstOfPairFlag(firstOfPair);
        read.setSecondOfPairFlag(!firstOfPair);

        FragmentCoords coords = FragmentCoords.fromRead(read, true);

        String actualCollapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);
        String expectedCollapsedKey = mateReversed ? "1:R:U" : "1:F:U";

        assertEquals(expectedCollapsedKey, actualCollapsedKey);
    }

    @Test
    public void testCollapseToKeyWithoutCoordinatesUnmappedMate()
    {
        checkCollapseToKeyWithoutCoordinatesUnmappedMate(true, true);
        checkCollapseToKeyWithoutCoordinatesUnmappedMate(true, false);
        checkCollapseToKeyWithoutCoordinatesUnmappedMate(false, true);
        checkCollapseToKeyWithoutCoordinatesUnmappedMate(false, false);
    }

    private void checkCollapseToKeyWithoutCoordinatesUnmappedMate(boolean firstOfPair, boolean readReversed)
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                "READ_001", CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION,
                readReversed, false, null, false, NO_CIGAR);
        read.setMateUnmappedFlag(true);
        read.setFirstOfPairFlag(firstOfPair);
        read.setSecondOfPairFlag(!firstOfPair);

        FragmentCoords coords = FragmentCoords.fromRead(read, true);

        String actualCollapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);
        String expectedCollapsedKey = readReversed ? "1:R" : "1:F";

        assertEquals(expectedCollapsedKey, actualCollapsedKey);
    }

    @Test
    public void testCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate()
    {
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(true, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(true, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(false, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(false, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(true, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(true, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(false, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(false, false, true);
    }

    private void checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(boolean firstOfPair, boolean primaryReversed, boolean supplementaryReversed)
    {
        char primaryOrientation = primaryReversed ? SUPP_NEG_STRAND : SUPP_POS_STRAND;
        SupplementaryReadData primaryAlignment = new SupplementaryReadData(CHR_1, 100, primaryOrientation, TEST_READ_CIGAR, 60);
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                "READ_001", CHR_2, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, supplementaryReversed,
                true, primaryAlignment, false, NO_CIGAR);
        read.setMateUnmappedFlag(true);
        read.setFirstOfPairFlag(firstOfPair);
        read.setSecondOfPairFlag(!firstOfPair);

        FragmentCoords coords = FragmentCoords.fromRead(read, true);

        String actualCollapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);
        String expectedCollapsedKey = primaryReversed ? "1:R:S" : "1:F:S";

        assertEquals(expectedCollapsedKey, actualCollapsedKey);
    }

    @Test
    public void testCollapseToKeyWithoutCoordinatesMappedPrimaryPair()
    {
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, true, true, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, true, false, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, true, true, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, true, false, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, true, true, false);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, true, false, false);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, true, true, false);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, true, false, false);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, false, true, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, false, false, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, false, true, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, false, false, true);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, false, true, false);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(true, false, false, false);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, false, true, false);
        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(false, false, false, false);
    }

    private void checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(boolean firstOfPair, boolean isLower, boolean readReversed, boolean mateReversed)
    {
        int readStart = isLower ? 100 : 1_000;
        int mateStart = isLower ? 1_000 : 100;
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                "READ_001", CHR_1, readStart, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, mateStart, readReversed, false, null,
                mateReversed, TEST_READ_CIGAR);
        read.setFirstOfPairFlag(firstOfPair);
        read.setSecondOfPairFlag(!firstOfPair);

        FragmentCoords coords = FragmentCoords.fromRead(read, true);

        String actualCollapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);

        String lowerOrientation = coords.OrientLower == FORWARD ? "F" : "R";
        String upperOrientation = coords.OrientUpper == FORWARD ? "F" : "R";
        String isLowerStr = isLower ? "L" : "U";
        String expectedCollapsedKey = format("1:%s:1:%s:%s", lowerOrientation, upperOrientation, isLowerStr);

        assertEquals(expectedCollapsedKey, actualCollapsedKey);
    }

    @Test
    public void testCollapseToKeyWithoutCoordinatesSupplementaryRead()
    {
        // TODO: should make variables to switch between values
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, true, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, false, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, true, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, false, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, true, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, false, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, true, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, false, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, true, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, false, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, true, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, false, true, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, true, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, false, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, true, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, false, true, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, true, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, false, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, true, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, false, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, true, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, true, false, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, true, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, true, false, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, true, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, false, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, true, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, false, false, true);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, true, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(true, false, false, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, true, false, false);
        checkCollapseToKeyWithoutCoordinatesSupplementaryRead(false, false, false, false, false);
    }

    private void checkCollapseToKeyWithoutCoordinatesSupplementaryRead(
            boolean firstOfPair, boolean isLower, boolean primaryReversed, boolean supplementaryReversed, boolean mateReversed)
    {
        int primaryStart = isLower ? 100 : 1_000;
        int mateStart = isLower ? 1_000 : 100;
        char primaryOrientation = primaryReversed ? SUPP_NEG_STRAND : SUPP_POS_STRAND;
        SupplementaryReadData primaryAlignment = new SupplementaryReadData(
                CHR_1, primaryStart, primaryOrientation, TEST_READ_CIGAR, 60);
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                "READ_001", CHR_2, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, mateStart, supplementaryReversed, true,
                primaryAlignment, mateReversed, TEST_READ_CIGAR);
        read.setFirstOfPairFlag(firstOfPair);
        read.setSecondOfPairFlag(!firstOfPair);

        FragmentCoords coords = FragmentCoords.fromRead(read, true);

        String actualCollapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);

        String lowerOrientation = coords.OrientLower == FORWARD ? "F" : "R";
        String upperOrientation = coords.OrientUpper == FORWARD ? "F" : "R";
        String isLowerStr = isLower ? "L" : "U";
        String expectedCollapsedKey = format("1:%s:1:%s:%s:S", lowerOrientation, upperOrientation, isLowerStr);

        assertEquals(expectedCollapsedKey, actualCollapsedKey);
    }

    @Test
    public void testCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos()
    {
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(true, false, false);
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(true, false, true);
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(true, true, false);
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(true, true, true);
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(false, false, false);
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(false, false, true);
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(false, true, false);
        checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(false, true, true);
    }

    private static final int FIVE_PRIME_POS = 200;

    private void checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(boolean firstOfPair, boolean readReversed, boolean mateReversed)
    {
        int readStart = !readReversed ? FIVE_PRIME_POS : FIVE_PRIME_POS - (TEST_READ_BASES.length() - 1);
        int mateStart = !mateReversed ? FIVE_PRIME_POS : FIVE_PRIME_POS - (TEST_READ_BASES.length() - 1);
        SAMRecord read = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, readStart, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                mateStart, readReversed, false, null, mateReversed, TEST_READ_CIGAR);
        read.setFirstOfPairFlag(firstOfPair);
        read.setSecondOfPairFlag(!firstOfPair);

        FragmentCoords coords = FragmentCoords.fromRead(read, true);

        assertEquals(FIVE_PRIME_POS, coords.PositionLower);
        assertEquals(FIVE_PRIME_POS, coords.PositionUpper);
        assertEquals(firstOfPair, coords.ReadIsLower);

        String actualCollapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);

        String lowerOrientation = coords.OrientLower == FORWARD ? "F" : "R";
        String upperOrientation = coords.OrientUpper == FORWARD ? "F" : "R";
        String isLowerStr = firstOfPair ? "L" : "U";
        String expectedCollapsedKey = format("1:%s:1:%s:%s", lowerOrientation, upperOrientation, isLowerStr);

        assertEquals(expectedCollapsedKey, actualCollapsedKey);
    }

}
