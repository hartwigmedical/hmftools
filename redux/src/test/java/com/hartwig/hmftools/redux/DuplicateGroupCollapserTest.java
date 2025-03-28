package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.collapseToKeyWithoutCoordinates;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.stream.Stream;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.common.FragmentCoords;

import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.TestFactory;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroupCollapserTest
{
    private static final String TEST_READ_BASES = "A".repeat(100);
    private static final String TEST_READ_CIGAR = "100M";

    @TestFactory
    public Stream<DynamicTest> testCollapseToKeyWithoutCoordinatesUnmappedRead()
    {
        record TestCase(String name, boolean firstOfPair, boolean mateReversed)
        {
            public void check()
            {
                SAMRecord read = SamRecordTestUtils.createSamRecord(
                        "READ_001", NO_CHROMOSOME_NAME, NO_POSITION, TEST_READ_BASES, NO_CIGAR, CHR_1, 1_000, false, false, null,
                        mateReversed, TEST_READ_CIGAR);
                read.setReadUnmappedFlag(true);
                read.setFirstOfPairFlag(firstOfPair);
                read.setSecondOfPairFlag(!firstOfPair);

                FragmentCoords coords = FragmentCoords.fromRead(read, true);

                String actualCollapsedKey = collapseToKeyWithoutCoordinates(coords, true);
                String expectedCollapsedKey = mateReversed ? "1:R:U" : "1:F:U";

                assertEquals(expectedCollapsedKey, actualCollapsedKey);
            }
        }

        TestCase[] testCases = new TestCase[] {
                new TestCase("firstOfPair_mateReverse", true, true),
                new TestCase("firstOfPair_mateForward", true, false),
                new TestCase("secondOfPair_mateReverse", false, true),
                new TestCase("secondOfPair_mateForward", false, false),
        };

        return DynamicTest.stream(Stream.of(testCases), TestCase::name, TestCase::check);
    }

    @TestFactory
    public Stream<DynamicTest> testCollapseToKeyWithoutCoordinatesUnmappedMate()
    {
        record TestCase(String name, boolean firstOfPair, boolean readReversed)
        {
            public void check()
            {
                SAMRecord read = SamRecordTestUtils.createSamRecord(
                        "READ_001", CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, readReversed, false,
                        null, false, NO_CIGAR);
                read.setMateUnmappedFlag(true);
                read.setFirstOfPairFlag(firstOfPair);
                read.setSecondOfPairFlag(!firstOfPair);

                FragmentCoords coords = FragmentCoords.fromRead(read, true);

                String actualCollapsedKey = collapseToKeyWithoutCoordinates(coords, true);
                String expectedCollapsedKey = readReversed ? "1:R" : "1:F";

                assertEquals(expectedCollapsedKey, actualCollapsedKey);
            }
        }

        TestCase[] testCases = new TestCase[] {
                new TestCase("firstOfPair_readReverse", true, true),
                new TestCase("firstOfPair_readForward", true, false),
                new TestCase("secondOfPair_readReverse", false, true),
                new TestCase("secondOfPair_readForward", false, false),
        };

        return DynamicTest.stream(Stream.of(testCases), TestCase::name, TestCase::check);
    }

    @TestFactory
    public Stream<DynamicTest> testCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate()
    {
        record TestCase(String name, boolean firstOfPair, boolean primaryReversed, boolean supplementaryReversed)
        {
            public void check()
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

                String actualCollapsedKey = collapseToKeyWithoutCoordinates(coords, true);
                String expectedCollapsedKey = primaryReversed ? "1:R:S" : "1:F:S";

                assertEquals(expectedCollapsedKey, actualCollapsedKey);
            }
        }

        TestCase[] testCases = new TestCase[] {
                new TestCase("firstOfPair_primaryReverse_suppForward", true, true, false),
                new TestCase("firstOfPair_primaryForward_suppForward", true, false, false),
                new TestCase("secondOfPair_primaryReverse_suppForward", false, true, false),
                new TestCase("secondOfPair_primaryForward_suppForward", false, false, false),
                new TestCase("firstOfPair_primaryReverse_suppReverse", true, true, true),
                new TestCase("firstOfPair_primaryForward_suppReverse", true, false, true),
                new TestCase("secondOfPair_primaryReverse_suppReverse", false, true, true),
                new TestCase("secondOfPair_primaryForward_suppReverse", false, false, true),
        };

        return DynamicTest.stream(Stream.of(testCases), TestCase::name, TestCase::check);
    }

    @TestFactory
    public Stream<DynamicTest> testCollapseToKeyWithoutCoordinatesMappedPrimaryPair()
    {
        record TestCase(String name, boolean firstOfPair, boolean isLower, boolean readReversed, boolean mateReversed)
        {
            public void check()
            {
                int readStart = isLower ? 100 : 1_000;
                int mateStart = isLower ? 1_000 : 100;
                SAMRecord read = SamRecordTestUtils.createSamRecord(
                        "READ_001", CHR_1, readStart, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, mateStart, readReversed, false, null,
                        mateReversed, TEST_READ_CIGAR);
                read.setFirstOfPairFlag(firstOfPair);
                read.setSecondOfPairFlag(!firstOfPair);

                FragmentCoords coords = FragmentCoords.fromRead(read, true);

                String actualCollapsedKey = collapseToKeyWithoutCoordinates(coords, true);

                String lowerOrientation = coords.OrientLower == FORWARD ? "F" : "R";
                String upperOrientation = coords.OrientUpper == FORWARD ? "F" : "R";
                String isLowerStr = isLower ? "L" : "U";
                String expectedCollapsedKey = format("1:%s:1:%s:%s", lowerOrientation, upperOrientation, isLowerStr);
                if(!coords.keyNonOriented().equals(coords.Key))
                    expectedCollapsedKey = expectedCollapsedKey + ":N";

                assertEquals(expectedCollapsedKey, actualCollapsedKey);
            }
        }

        TestCase[] testCases = new TestCase[] {
                new TestCase("firstOfPair_isLower_readReverse_mateReverse", true, true, true, true),
                new TestCase("firstOfPair_isLower_readForward_mateReverse", true, true, false, true),
                new TestCase("secondOfPair_isLower_readReverse_mateReverse", false, true, true, true),
                new TestCase("secondOfPair_isLower_readForward_mateReverse", false, true, false, true),
                new TestCase("firstOfPair_isLower_readReverse_mateForward", true, true, true, false),
                new TestCase("firstOfPair_isLower_readForward_mateForward", true, true, false, false),
                new TestCase("secondOfPair_isLower_readReverse_mateForward", false, true, true, false),
                new TestCase("secondOfPair_isLower_readForward_mateForward", false, true, false, false),
                new TestCase("firstOfPair_isUpper_readReverse_mateReverse", true, false, true, true),
                new TestCase("firstOfPair_isUpper_readForward_mateReverse", true, false, false, true),
                new TestCase("secondOfPair_isUpper_readReverse_mateReverse", false, false, true, true),
                new TestCase("secondOfPair_isUpper_readForward_mateReverse", false, false, false, true),
                new TestCase("firstOfPair_isUpper_readReverse_mateForward", true, false, true, false),
                new TestCase("firstOfPair_isUpper_readForward_mateForward", true, false, false, false),
                new TestCase("secondOfPair_isUpper_readReverse_mateForward", false, false, true, false),
                new TestCase("secondOfPair_isUpper_readForward_mateForward", false, false, false, false),
        };

        return DynamicTest.stream(Stream.of(testCases), TestCase::name, TestCase::check);
    }

    @TestFactory
    public Stream<DynamicTest> testCollapseToKeyWithoutCoordinatesSupplementaryRead()
    {
        record TestCase(String name, boolean firstOfPair, boolean isLower, boolean primaryReversed, boolean supplementaryReversed,
                        boolean mateReversed)
        {
            public void check()
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

                String actualCollapsedKey = collapseToKeyWithoutCoordinates(coords, true);

                String lowerOrientation = coords.OrientLower == FORWARD ? "F" : "R";
                String upperOrientation = coords.OrientUpper == FORWARD ? "F" : "R";
                String isLowerStr = isLower ? "L" : "U";
                String expectedCollapsedKey = format("1:%s:1:%s:%s:S", lowerOrientation, upperOrientation, isLowerStr);
                if(!coords.keyNonOriented().equals(coords.Key))
                    expectedCollapsedKey = expectedCollapsedKey + ":N";

                assertEquals(expectedCollapsedKey, actualCollapsedKey);
            }
        }

        TestCase[] testCases = new TestCase[] {
                new TestCase("firstOfPair_isLower_primaryReverse_supplementaryReverse_mateReverse", true, true, true, true, true),
                new TestCase("firstOfPair_isLower_primaryForward_supplementaryReverse_mateReverse", true, true, false, true, true),
                new TestCase("secondOfPair_isLower_primaryReverse_supplementaryReverse_mateReverse", false, true, true, true, true),
                new TestCase("secondOfPair_isLower_primaryForward_supplementaryReverse_mateReverse", false, true, false, true, true),
                new TestCase("firstOfPair_isLower_primaryReverse_supplementaryReverse_mateForward", true, true, true, true, false),
                new TestCase("firstOfPair_isLower_primaryForward_supplementaryReverse_mateForward", true, true, false, true, false),
                new TestCase("secondOfPair_isLower_primaryReverse_supplementaryReverse_mateForward", false, true, true, true, false),
                new TestCase("secondOfPair_isLower_primaryForward_supplementaryReverse_mateForward", false, true, false, true, false),
                new TestCase("firstOfPair_isUpper_primaryReverse_supplementaryReverse_mateReverse", true, false, true, true, true),
                new TestCase("firstOfPair_isUpper_primaryForward_supplementaryReverse_mateReverse", true, false, false, true, true),
                new TestCase("secondOfPair_isUpper_primaryReverse_supplementaryReverse_mateReverse", false, false, true, true, true),
                new TestCase("secondOfPair_isUpper_primaryForward_supplementaryReverse_mateReverse", false, false, false, true, true),
                new TestCase("firstOfPair_isUpper_primaryReverse_supplementaryReverse_mateForward", true, false, true, true, false),
                new TestCase("firstOfPair_isUpper_primaryForward_supplementaryReverse_mateForward", true, false, false, true, false),
                new TestCase("secondOfPair_isUpper_primaryReverse_supplementaryReverse_mateForward", false, false, true, true, false),
                new TestCase("secondOfPair_isUpper_primaryForward_supplementaryReverse_mateForward", false, false, false, true, false),
                new TestCase("firstOfPair_isLower_primaryReverse_supplementaryForward_mateReverse", true, true, true, false, true),
                new TestCase("firstOfPair_isLower_primaryForward_supplementaryForward_mateReverse", true, true, false, false, true),
                new TestCase("secondOfPair_isLower_primaryReverse_supplementaryForward_mateReverse", false, true, true, false, true),
                new TestCase("secondOfPair_isLower_primaryForward_supplementaryForward_mateReverse", false, true, false, false, true),
                new TestCase("firstOfPair_isLower_primaryReverse_supplementaryForward_mateForward", true, true, true, false, false),
                new TestCase("firstOfPair_isLower_primaryForward_supplementaryForward_mateForward", true, true, false, false, false),
                new TestCase("secondOfPair_isLower_primaryReverse_supplementaryForward_mateForward", false, true, true, false, false),
                new TestCase("secondOfPair_isLower_primaryForward_supplementaryForward_mateForward", false, true, false, false, false),
                new TestCase("firstOfPair_isUpper_primaryReverse_supplementaryForward_mateReverse", true, false, true, false, true),
                new TestCase("firstOfPair_isUpper_primaryForward_supplementaryForward_mateReverse", true, false, false, false, true),
                new TestCase("secondOfPair_isUpper_primaryReverse_supplementaryForward_mateReverse", false, false, true, false, true),
                new TestCase("secondOfPair_isUpper_primaryForward_supplementaryForward_mateReverse", false, false, false, false, true),
                new TestCase("firstOfPair_isUpper_primaryReverse_supplementaryForward_mateForward", true, false, true, false, false),
                new TestCase("firstOfPair_isUpper_primaryForward_supplementaryForward_mateForward", true, false, false, false, false),
                new TestCase("secondOfPair_isUpper_primaryReverse_supplementaryForward_mateForward", false, false, true, false, false),
                new TestCase("secondOfPair_isUpper_primaryForward_supplementaryForward_mateForward", false, false, false, false, false),
        };

        return DynamicTest.stream(Stream.of(testCases), TestCase::name, TestCase::check);
    }
}
