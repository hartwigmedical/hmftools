package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.common.DuplicateGroupBuilder.findDuplicateFragments;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.redux.common.FragmentUtils.calcFragmentStatus;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createFragment;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.redux.common.CandidateDuplicates;
import com.hartwig.hmftools.redux.common.Fragment;

import org.junit.Test;

public class FragmentClassificationTest
{
    private final ReadIdGenerator mReadIdGen;

    public FragmentClassificationTest()
    {
        mReadIdGen = new ReadIdGenerator();
    }

    private void initialiseFragmentCoordinates(final Fragment fragment) { fragment.intialiseCoordinates(true); }

    @Test
    public void testFragmentPairStatus()
    {
        Fragment frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);
        initialiseFragmentCoordinates(frag1);

        Fragment frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, false, null);
        initialiseFragmentCoordinates(frag2);

        assertEquals(CANDIDATE, calcFragmentStatus(frag1, frag2, true));

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                true, false, null);
        initialiseFragmentCoordinates(frag2);

        assertEquals(NONE, calcFragmentStatus(frag1, frag2, true));

        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                true, false, null);
        initialiseFragmentCoordinates(frag1);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                true, false, null);
        initialiseFragmentCoordinates(frag2);

        assertEquals(CANDIDATE, calcFragmentStatus(frag1, frag2, true));

        // diff positions at end
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                true, false, null);
        initialiseFragmentCoordinates(frag1);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, "101M", CHR_1, 201,
                true, false, null);
        initialiseFragmentCoordinates(frag2);

        assertEquals(NONE, calcFragmentStatus(frag1, frag2, true));

        // diff mate orientations
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false,
                CHR_1, 200, false, TEST_READ_CIGAR);
        initialiseFragmentCoordinates(frag1);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_CIGAR, false,
                CHR_1, 200, true, TEST_READ_CIGAR);
        initialiseFragmentCoordinates(frag2);

        assertEquals(NONE, calcFragmentStatus(frag1, frag2, true));

        // unpaired matching
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, "", 0,
                false, false, null);
        initialiseFragmentCoordinates(frag1);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, "", 0,
                false, false, null);
        initialiseFragmentCoordinates(frag2);

        assertEquals(DUPLICATE, calcFragmentStatus(frag1, frag2, true));

        frag2.reads().get(0).setInferredInsertSize(100);
        initialiseFragmentCoordinates(frag2);

        assertEquals(NONE, calcFragmentStatus(frag1, frag2, true));

        // mates present and matching
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1000,
                false, false, null);
        initialiseFragmentCoordinates(frag1);

        frag1.addRead(createSamRecord(mReadIdGen.currentId(), CHR_1, 1000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null));

        assertTrue(frag1.primaryReadsPresent());

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, false, null);
        initialiseFragmentCoordinates(frag2);

        frag2.addRead(createSamRecord(mReadIdGen.currentId(), CHR_1, 1000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null));

        assertTrue(frag2.primaryReadsPresent());

        assertEquals(DUPLICATE, calcFragmentStatus(frag1, frag2, true));
    }

    @Test
    public void testClassifyPositionFragment()
    {
        // test 1: single group of candidates with one NONE

        Fragment frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);

        Fragment frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, false, null);

        // different chromosome
        Fragment frag3 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 201,
                false, false, null);

        List<Fragment> fragments = Lists.newArrayList(frag1, frag2, frag3);

        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList();
        List<List<Fragment>> duplicateGroups = Lists.newArrayList();

        findDuplicateFragments(fragments, resolvedFragments, duplicateGroups, candidateDuplicatesList, false);

        assertEquals(1, resolvedFragments.size());
        assertEquals(frag3, resolvedFragments.get(0));

        assertEquals(1, candidateDuplicatesList.size());
        CandidateDuplicates candidateDuplicates = candidateDuplicatesList.get(0);
        assertEquals(2, candidateDuplicates.fragmentCount());
        assertTrue(candidateDuplicates.fragments().contains(frag1));
        assertTrue(candidateDuplicates.fragments().contains(frag2));

        // order doesn't matter
        fragments = Lists.newArrayList(frag3, frag1, frag2);

        resolvedFragments.clear();
        candidateDuplicatesList.clear();
        duplicateGroups.clear();

        findDuplicateFragments(fragments, resolvedFragments, duplicateGroups, candidateDuplicatesList, false);

        assertEquals(1, resolvedFragments.size());
        assertEquals(frag3, resolvedFragments.get(0));

        assertEquals(1, candidateDuplicatesList.size());
        candidateDuplicates = candidateDuplicatesList.get(0);
        assertEquals(2, candidateDuplicates.fragmentCount());
        assertTrue(candidateDuplicates.fragments().contains(frag1));
        assertTrue(candidateDuplicates.fragments().contains(frag2));

        mReadIdGen.reset();

        // 2 separate candidate groups, and one which requires a link between latter fragments
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 260,
                false, false, null);

        frag3 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 310,
                false, false, null);

        Fragment frag4 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 400,
                false, false, null);

        Fragment frag5 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, 400,
                false, false, null);

        fragments = Lists.newArrayList(frag4, frag1, frag2, frag5, frag3);

        resolvedFragments.clear();
        candidateDuplicatesList.clear();

        findDuplicateFragments(fragments, resolvedFragments, duplicateGroups, candidateDuplicatesList, false);

        assertEquals(0, resolvedFragments.size());

        assertEquals(2, candidateDuplicatesList.size());
        candidateDuplicates = candidateDuplicatesList.stream().filter(x -> x.fragmentCount() == 2).findFirst().orElse(null);
        assertNotNull(candidateDuplicates);
        assertTrue(candidateDuplicates.fragments().contains(frag4));
        assertTrue(candidateDuplicates.fragments().contains(frag5));

        candidateDuplicates = candidateDuplicatesList.stream().filter(x -> x.fragmentCount() == 3).findFirst().orElse(null);
        assertNotNull(candidateDuplicates);        assertTrue(candidateDuplicates.fragments().contains(frag1));
        assertTrue(candidateDuplicates.fragments().contains(frag2));
        assertTrue(candidateDuplicates.fragments().contains(frag3));

        mReadIdGen.reset();

        // more complicate example of successive links
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);

        frag3 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 260,
                false, false, null);

        frag4 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 260,
                false, false, null);

        frag5 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 310,
                false, false, null);

        Fragment frag6 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 400,
                false, false, null);

        Fragment frag7 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 400,
                false, false, null);

        fragments = Lists.newArrayList(frag1, frag2, frag6, frag7, frag3, frag4, frag5);

        resolvedFragments.clear();
        candidateDuplicatesList.clear();
        duplicateGroups.clear();

        findDuplicateFragments(fragments, resolvedFragments, duplicateGroups, candidateDuplicatesList, false);

        assertEquals(0, resolvedFragments.size());
        assertEquals(1, candidateDuplicatesList.size());
        candidateDuplicates = candidateDuplicatesList.get(0);
        assertEquals(7, candidateDuplicates.fragmentCount());
        assertTrue(candidateDuplicates.fragments().stream().allMatch(x -> x.status() == CANDIDATE));

        // and again, based on an observed example from RNA
        mReadIdGen.reset();

        // more complicate example of successive links
        frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 8108,
                false, false, null);

        frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 7389,
                false, false, null);

        frag3 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 7274,
                false, false, null);

        frag4 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 8117,
                false, false, null);

        frag5 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 7338,
                false, false, null);

        frag6 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 8117,
                false, false, null);

        frag7 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 7262,
                false, false, null);

        Fragment frag8 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 7261,
                false, false, null);

        Fragment frag9 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 4295,
                false, false, null);

        fragments = Lists.newArrayList(frag1, frag2, frag3, frag4, frag5, frag6, frag7, frag8, frag9);

        resolvedFragments.clear();
        candidateDuplicatesList.clear();
        duplicateGroups.clear();

        findDuplicateFragments(fragments, resolvedFragments, duplicateGroups, candidateDuplicatesList, false);

        assertEquals(1, resolvedFragments.size());
        assertTrue(resolvedFragments.contains(frag9));

        assertEquals(2, candidateDuplicatesList.size());
        candidateDuplicates = candidateDuplicatesList.stream().filter(x -> x.fragmentCount() == 3).findFirst().orElse(null);
        assertNotNull(candidateDuplicates);
        assertTrue(candidateDuplicates.fragments().contains(frag1));
        assertTrue(candidateDuplicates.fragments().contains(frag4));
        assertTrue(candidateDuplicates.fragments().contains(frag6));
        assertTrue(candidateDuplicates.fragments().stream().allMatch(x -> x.status() == CANDIDATE));

        candidateDuplicates = candidateDuplicatesList.stream().filter(x -> x.fragmentCount() == 5).findFirst().orElse(null);
        assertNotNull(candidateDuplicates);

        assertTrue(candidateDuplicates.fragments().contains(frag2));
        assertTrue(candidateDuplicates.fragments().contains(frag3));
        assertTrue(candidateDuplicates.fragments().contains(frag5));
        assertTrue(candidateDuplicates.fragments().contains(frag7));
        assertTrue(candidateDuplicates.fragments().contains(frag8));
    }
}
