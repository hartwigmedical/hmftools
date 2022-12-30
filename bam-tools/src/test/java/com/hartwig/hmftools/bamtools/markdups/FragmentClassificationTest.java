package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.checkDuplicateFragments;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragment;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.ReadIdGenerator;

import org.junit.Test;

public class FragmentClassificationTest
{
    private final ReadIdGenerator mReadIdGen;

    public FragmentClassificationTest()
    {
        mReadIdGen = new ReadIdGenerator();
    }

    /*
    @Test
    public void testFragmentPairStatus()
    {
        Fragment frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);

        Fragment frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, false, null);

        assertEquals(UNCLEAR, calcFragmentStatus(frag1, frag2));

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                true, false, null);

        assertEquals(NONE, calcFragmentStatus(frag1, frag2));

        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                true, false, null);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                true, false, null);

        assertEquals(UNCLEAR, calcFragmentStatus(frag1, frag2));

        // diff positions at end
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                true, false, null);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "101M", CHR_1, 201,
                true, false, null);

        assertEquals(NONE, calcFragmentStatus(frag1, frag2));

        // unpaired matching
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, "", 0,
                false, false, null);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, "", 0,
                false, false, null);

        assertEquals(DUPLICATE, calcFragmentStatus(frag1, frag2));

        // mates present and matching
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1000,
                false, false, null);

        frag1.addRead(createSamRecord(mReadIdGen.currentId(), CHR_1, 1000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null));

        assertTrue(frag1.primaryReadsPresent());

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, false, null);

        frag2.addRead(createSamRecord(mReadIdGen.currentId(), CHR_1, 1000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null));

        assertTrue(frag2.primaryReadsPresent());

        assertEquals(DUPLICATE, calcFragmentStatus(frag1, frag2));
    }
    */

    @Test
    public void testFragmentDuplicates()
    {
        List<Fragment> fragments = Lists.newArrayList();

        Fragment read1 = createFragment(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false, CHR_1, 300, true, TEST_READ_CIGAR);

        Fragment dupFragment = checkDuplicateFragments(read1, fragments);
        assertNull(dupFragment);
        assertEquals(1, fragments.size());

        // different from the first
        Fragment read2 = createFragment(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false, CHR_1, 400, true, TEST_READ_CIGAR);
        setBaseQualities(read2, DEFAULT_QUAL - 1);

        dupFragment = checkDuplicateFragments(read2, fragments);
        assertNull(dupFragment);
        assertEquals(2, fragments.size());

        // duplicate of the first, lower qual
        Fragment read3 = createFragment(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false, CHR_1, 300, true, TEST_READ_CIGAR);
        setBaseQualities(read3, DEFAULT_QUAL - 1);

        dupFragment = checkDuplicateFragments(read3, fragments);
        assertNotNull(dupFragment);
        assertEquals(2, fragments.size());
        assertEquals(DUPLICATE, dupFragment.status());

        // duplicate of the second, higher qual so swaps out
        Fragment read4 = createFragment(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false, CHR_1, 400, true, TEST_READ_CIGAR);

        dupFragment = checkDuplicateFragments(read4, fragments);
        assertEquals(dupFragment, read2);
        assertEquals(2, fragments.size());

        // no duplicates
        Fragment read5 = createFragment(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false, CHR_2, 300, true, TEST_READ_CIGAR);

        dupFragment = checkDuplicateFragments(read5, fragments);
        assertNull(dupFragment);

        assertEquals(PRIMARY, read1.status());
        assertEquals(PRIMARY, read4.status());
        assertEquals(UNSET, read5.status());
        assertEquals(1, read1.duplicateCount());
        assertEquals(1, read4.duplicateCount());

        // duplicates of the second, lower qual
        Fragment read6 = createFragment(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false, CHR_1, 400, true, TEST_READ_CIGAR);
        setBaseQualities(read6, DEFAULT_QUAL - 2);

        dupFragment = checkDuplicateFragments(read6, fragments);
        assertEquals(dupFragment, read6);

        Fragment read7 = createFragment(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false, CHR_1, 400, true, TEST_READ_CIGAR);
        setBaseQualities(read7, DEFAULT_QUAL - 2);

        dupFragment = checkDuplicateFragments(read7, fragments);
        assertEquals(dupFragment, read7);

        assertEquals(PRIMARY, read4.status());
        assertEquals(3, read4.duplicateCount());
    }
}
