package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.redux.old.FragmentUtils.calcFragmentStatus;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createFragment;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.redux.old.FragmentOld;

import org.junit.Test;

public class FragClassificationOldTest
{
    private final ReadIdGenerator mReadIdGen;

    public FragClassificationOldTest()
    {
        mReadIdGen = new ReadIdGenerator();
    }

    private void initialiseFragmentCoordinates(final FragmentOld fragment) { fragment.intialiseCoordinates(true); }

    @Test
    public void testFragmentPairStatus()
    {
        FragmentOld frag1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);
        initialiseFragmentCoordinates(frag1);

        FragmentOld frag2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
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
}
