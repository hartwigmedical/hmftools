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
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createFragmentCoords;
import static com.hartwig.hmftools.redux.TestUtils.createTestConfig;
import static com.hartwig.hmftools.redux.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.redux.duplicate.CollapseUtils.collapseToNonOrientedKeyWithoutCoordinates;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_FORWARD_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_REVERSE_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_LOWER_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_UPPER_STR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroupTest
{
    @Test
    public void testNonConsensusPrimaryRead()
    {
        ReduxConfig config = createTestConfig();

        DuplicateGroupBuilder duplicateGroupBuilder = new DuplicateGroupBuilder(config);

        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();

        int readPos = 100;
        int matePos = 200;

        SAMRecord read1 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, matePos, false, false, null, true, TEST_READ_CIGAR);

        setBaseQualities(read1, 25);

        FragmentCoords fragmentCoords = createFragmentCoords(read1);

        SAMRecord read2 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, matePos, false, false, null, true, TEST_READ_CIGAR);

        setBaseQualities(read2, 25);

        SAMRecord read3 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, matePos, false, false, null, true, TEST_READ_CIGAR);

        setBaseQualities(read3, 11);

        DuplicateGroup duplicateGroup = new DuplicateGroup(List.of(read1, read2, read3), fragmentCoords);
        duplicateGroups.add(duplicateGroup);

        duplicateGroupBuilder.processDuplicateGroups(duplicateGroups, Collections.emptyList(), true);

        assertTrue(duplicateGroup.isPrimaryRead(read1));
    }

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
        String expectedCollapsedKey = mateReversed ? "1:R:N" : "1:F:N";

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
        for(int a = 0; a <= 1; ++a)
        {
            boolean firstInPair = (a == 0);

            for(int b = 0; b <= 1; ++b)
            {
                boolean primaryReversed = (b == 0);

                for(int c = 0; c <= 1; ++c)
                {
                    boolean suppReversed = (c == 0);

                    checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(firstInPair, primaryReversed, suppReversed);
                }
            }
        }
    }

    private void checkCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate(
            boolean firstOfPair, boolean primaryReversed, boolean supplementaryReversed)
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
        for(int a = 0; a <= 1; ++a)
        {
            boolean firstInPair = (a == 0);

            for(int b = 0; b <= 1; ++b)
            {
                boolean readIsLower = (b == 0);

                for(int c = 0; c <= 1; ++c)
                {
                    boolean readReversed = (c == 0);

                    for(int d = 0; d <= 1; ++d)
                    {
                        boolean mateReversed = (d == 0);
                        checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(firstInPair, readIsLower, readReversed, mateReversed);
                    }
                }
            }
        }
    }

    private void checkCollapseToKeyWithoutCoordinatesMappedPrimaryPair(
            boolean firstOfPair, boolean isLower, boolean readReversed, boolean mateReversed)
    {
        int readStart = isLower ? 100 : 1_000;
        int mateStart = isLower ? 1_000 : 100;
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                "READ_001", CHR_1, readStart, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, mateStart, readReversed,
                false, null, mateReversed, TEST_READ_CIGAR);
        read.setFirstOfPairFlag(firstOfPair);
        read.setSecondOfPairFlag(!firstOfPair);

        FragmentCoords coords = FragmentCoords.fromRead(read, true);

        String actualCollapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);

        String lowerOrientation = coords.OrientLower == FORWARD ? COORD_ORIENT_FORWARD_STR : COORD_ORIENT_REVERSE_STR;
        String upperOrientation = coords.OrientUpper == FORWARD ? COORD_ORIENT_FORWARD_STR : COORD_ORIENT_REVERSE_STR;
        String isLowerStr = isLower ? COORD_READ_LOWER_STR : COORD_READ_UPPER_STR;
        String expectedCollapsedKey = format("1:%s:1:%s:%s", lowerOrientation, upperOrientation, isLowerStr);

        assertEquals(expectedCollapsedKey, actualCollapsedKey);
    }

    @Test
    public void testCollapseToKeyWithoutCoordinatesSupplementaryRead()
    {
        for(int a = 0; a <= 1; ++a)
        {
            boolean firstInPair = (a == 0);

            for(int b = 0; b <= 1; ++b)
            {
                boolean readIsLower = (b == 0);

                for(int c = 0; c <= 1; ++c)
                {
                    boolean primaryReversed = (c == 0);

                    for(int d = 0; d <= 1; ++d)
                    {
                        boolean suppReversed = (d == 0);

                        for(int e = 0; e <= 1; ++e)
                        {
                            boolean mateReversed = (e == 0);
                            checkCollapseToKeyWithoutCoordinatesSupplementaryRead(
                                    firstInPair, readIsLower, primaryReversed, suppReversed, mateReversed);
                        }
                    }
                }
            }
        }
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
        for(int a = 0; a <= 1; ++a)
        {
            boolean firstInPair = (a == 0);

            for(int b = 0; b <= 1; ++b)
            {
                boolean readReversed = (b == 0);

                for(int c = 0; c <= 1; ++c)
                {
                    boolean mateReversed = (c == 0);
                    checkCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos(firstInPair, readReversed, mateReversed);
                }
            }
        }
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
