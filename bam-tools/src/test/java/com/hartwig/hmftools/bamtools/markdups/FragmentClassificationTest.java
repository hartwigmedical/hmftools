package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.calcFragmentStatus;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragmentsOld;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragment;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragmentPair;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createSamRecord;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.setBaseQualities;
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

        positionFragmentsList.add(frag1);
        positionFragmentsList.add(frag3);
        positionFragmentsList.add(frag2);

        classifyFragments(positionFragmentsList, resolvedFragments, incompletePositionFragments);

        assertEquals(1, resolvedFragments.size());
        assertEquals(1, incompletePositionFragments.size());

        assertTrue(resolvedFragments.contains(frag3));

        PositionFragments positionFragments = incompletePositionFragments.get(0);
        assertNotNull(positionFragments);
        assertTrue(positionFragments.Fragments.contains(frag1));
        assertTrue(positionFragments.Fragments.contains(frag2));

        assertEquals(UNCLEAR, frag1.status());
        assertEquals(UNCLEAR, frag2.status());
        assertEquals(NONE, frag3.status());

        // on the reverse strand
        mReadIdGen.reset();

        Fragment frag4 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 200,
                false, false, null); // unmatched since mate is elsewhere

        // matched on reverse strand
        Fragment frag5 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "100M", CHR_3, 200,
                true, false, null);

        Fragment frag6 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "95M5S", CHR_3, 199,
                true, false, null); // unmatched since mate is elsewhere

        positionFragmentsList = Lists.newArrayList(frag6, frag5, frag4);
        resolvedFragments.clear();
        incompletePositionFragments.clear();

        classifyFragments(positionFragmentsList, resolvedFragments, incompletePositionFragments);

        assertEquals(1, resolvedFragments.size());
        assertEquals(1, incompletePositionFragments.size());

        assertTrue(resolvedFragments.contains(frag4));

        positionFragments = incompletePositionFragments.get(0);
        assertNotNull(positionFragments);
        assertTrue(positionFragments.Fragments.contains(frag5));
        assertTrue(positionFragments.Fragments.contains(frag6));

        assertEquals(NONE, frag4.status());
        assertEquals(UNCLEAR, frag5.status());
        assertEquals(UNCLEAR, frag6.status());

        // now with some duplicates
        mReadIdGen.reset();

        // the first 3 are duplicates
        frag1 = createFragmentPair(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200, false);
        setBaseQualities(frag1, DEFAULT_QUAL - 2);

        frag2 = createFragmentPair(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200, false);
        setBaseQualities(frag2, DEFAULT_QUAL - 3);

        frag3 = createFragmentPair(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,false);
        setBaseQualities(frag5, DEFAULT_QUAL - 1);

        // then un-related
        frag4 = createFragmentPair(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 2000,false); // unmatched since mate is elsewhere

        // then 2 more duplicates
        frag5 = createFragmentPair(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, 200, false);
        setBaseQualities(frag5, DEFAULT_QUAL - 1);

        frag6 = createFragmentPair(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, 200, false); // unmatched since mate is elsewhere

        positionFragmentsList = Lists.newArrayList(frag1, frag5, frag4, frag6, frag3, frag2);
        resolvedFragments.clear();
        incompletePositionFragments.clear();

        classifyFragments(positionFragmentsList, resolvedFragments, incompletePositionFragments);

        assertEquals(6, resolvedFragments.size());

        assertEquals(DUPLICATE, frag1.status());
        assertEquals(DUPLICATE, frag2.status());
        assertEquals(DUPLICATE, frag3.status());
        assertEquals(NONE, frag4.status());
        assertEquals(DUPLICATE, frag5.status());
        assertEquals(PRIMARY, frag6.status());
    }
}
