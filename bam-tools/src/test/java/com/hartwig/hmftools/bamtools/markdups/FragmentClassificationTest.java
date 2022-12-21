package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.calcFragmentStatus;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragment;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

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

    @Test
    public void testFragmentClassification()
    {
        List<Fragment> positionFragmentsList = Lists.newArrayList();
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<PositionFragments> incompletePositionFragments = Lists.newArrayList();

        // a single fragment must be resolved as NONE
        Fragment frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100);
        positionFragmentsList.add(frag1);

        classifyFragments(positionFragmentsList, resolvedFragments, incompletePositionFragments);

        assertEquals(1, resolvedFragments.size());
        assertEquals(NONE, frag1.status());

        positionFragmentsList.clear();
        resolvedFragments.clear();
        incompletePositionFragments.clear();

        // a collection of fragments without their mates, all NONE or UNCLEAR
        mReadIdGen.reset();

        // paired without mates
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
            false, false, null);

        Fragment frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, false, null);

        // unmatched since mate is elsewhere
        Fragment frag3 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 2000,
                false, false, null);

        Fragment frag4 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 200,
                false, false, null); // unmatched since mate is elsewhere

        // matched on reverse strand
        Fragment frag5 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "100M", CHR_3, 200,
                true, false, null);

        Fragment frag6 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "95M5S", CHR_3, 199,
                true, false, null); // unmatched since mate is elsewhere

        positionFragmentsList.add(frag1);
        positionFragmentsList.add(frag6);
        positionFragmentsList.add(frag3);
        positionFragmentsList.add(frag5);
        positionFragmentsList.add(frag2);
        positionFragmentsList.add(frag4);

        classifyFragments(positionFragmentsList, resolvedFragments, incompletePositionFragments);

        assertEquals(2, resolvedFragments.size());
        assertEquals(2, incompletePositionFragments.size());

        assertTrue(resolvedFragments.contains(frag3));
        assertTrue(resolvedFragments.contains(frag4));

        PositionFragments positionFragments = incompletePositionFragments.stream().filter(x -> x.Position == 100).findFirst().orElse(null);
        assertNotNull(positionFragments);
        assertTrue(positionFragments.Fragments.contains(frag1));
        assertTrue(positionFragments.Fragments.contains(frag2));

        positionFragments = incompletePositionFragments.stream().filter(x -> x.Position == -199).findFirst().orElse(null);
        assertNotNull(positionFragments);
        assertTrue(positionFragments.Fragments.contains(frag5));
        assertTrue(positionFragments.Fragments.contains(frag6));

        assertEquals(UNCLEAR, frag1.status());
        assertEquals(UNCLEAR, frag2.status());
        assertEquals(NONE, frag3.status());
        assertEquals(NONE, frag4.status());
        assertEquals(UNCLEAR, frag5.status());
        assertEquals(UNCLEAR, frag6.status());

        positionFragmentsList.clear();
        resolvedFragments.clear();
        incompletePositionFragments.clear();
    }


}
