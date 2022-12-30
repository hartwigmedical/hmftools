package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.calcBaseQualAverage;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.findPrimaryFragment;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.getFragmentCoordinates;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.getUnclippedPosition;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragment;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createSamRecord;

import com.google.common.collect.Lists;

import static com.hartwig.hmftools.bamtools.markdups.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class FragmentUtilsTest
{
    @Test
    public void testFragmentCoords()
    {
        SAMRecord read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                false, false, null);

        int ucPos = getUnclippedPosition(read);
        assertEquals(100, ucPos);

        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "5S95M", CHR_1, 200,
                false, false, null);

        ucPos = getUnclippedPosition(read);
        assertEquals(95, ucPos);

        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "5S80M15S", CHR_1, 200,
                true, false, null);

        ucPos = getUnclippedPosition(read);
        assertEquals(194, ucPos);

        String cigarStr = "5S100M5S";
        ucPos = getUnclippedPosition(100, cigarStr, true);
        assertEquals(95, ucPos);

        ucPos = getUnclippedPosition(100, cigarStr, false);
        assertEquals(204, ucPos);

        cigarStr = "5S10M5D10I15D20M5S"; // +10 + 5, ignore 10, +15, +20, +5 = 40
        ucPos = getUnclippedPosition(100, cigarStr, false);
        assertEquals(154, ucPos);

        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                false, false, null);
        read.setMateNegativeStrandFlag(true);
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        FragmentCoordinates fragmentCoords = getFragmentCoordinates(read, true);
        assertEquals("1_100_1_299_R", fragmentCoords.Key);
        assertEquals(100, fragmentCoords.InitialPosition);

        // mate on earlier chromosome
        read = createSamRecord(TEST_READ_ID, CHR_2, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                false, false, null);
        read.setMateNegativeStrandFlag(false);
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        fragmentCoords = getFragmentCoordinates(read, true);
        assertEquals("1_200_2_100", fragmentCoords.Key);
        assertEquals(200, fragmentCoords.InitialPosition);

        fragmentCoords = getFragmentCoordinates(read, false);
        assertEquals("2_100_1_200", fragmentCoords.Key);
        assertEquals(100, fragmentCoords.InitialPosition);

        // mate in earlier position
        read = createSamRecord(TEST_READ_ID, CHR_1, 200, TEST_READ_BASES, "100M", CHR_1, 100,
                false, false, null);
        read.setMateNegativeStrandFlag(true);
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        fragmentCoords = getFragmentCoordinates(read, true);
        assertEquals("1_199_R_1_200", fragmentCoords.Key);
        assertEquals(-199, fragmentCoords.InitialPosition);
    }

    @Test
    public void testPrimaryDuplicateIdentification()
    {
        Fragment fragment = createFragment(TEST_READ_ID, CHR_1, 100);
        double baseAvg = calcBaseQualAverage(fragment);
        assertEquals(DEFAULT_QUAL, baseAvg, 0.1);

        Fragment fragment2 = createFragment(TEST_READ_ID, CHR_1, 100);

        SAMRecord read2 = fragment2.reads().get(0);
        setBaseQualities(read2, DEFAULT_QUAL - 1);

        assertEquals(DEFAULT_QUAL - 1, calcBaseQualAverage(fragment2), 0.1);

        List<Fragment> fragments = Lists.newArrayList(fragment2, fragment);

        Fragment primary = findPrimaryFragment(fragments, false);
        assertEquals(primary, fragment);
    }
}
