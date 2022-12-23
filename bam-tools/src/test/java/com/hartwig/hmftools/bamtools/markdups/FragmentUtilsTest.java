package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.calcBaseQualTotal;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.findPrimaryFragment;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.getUnclippedPosition;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragment;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createSamRecord;

import com.google.common.collect.Lists;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

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
    }

    @Test
    public void testPrimaryDuplicateIdentification()
    {
        Fragment fragment = createFragment(TEST_READ_ID, CHR_1, 100);
        int baseAvg = calcBaseQualTotal(fragment);
        assertEquals(DEFAULT_QUAL, baseAvg);

        Fragment fragment2 = createFragment(TEST_READ_ID, CHR_1, 100);

        SAMRecord read2 = fragment2.reads().get(0);
        for(int i = 0; i < read2.getReadLength(); ++i)
        {
            read2.getBaseQualities()[i] = DEFAULT_QUAL - 1;
        }

        assertEquals(DEFAULT_QUAL - 1, calcBaseQualTotal(fragment2));

        List<Fragment> fragments = Lists.newArrayList(fragment2, fragment);

        Fragment primary = findPrimaryFragment(fragments, false);
        assertEquals(primary, fragment);
    }
}
