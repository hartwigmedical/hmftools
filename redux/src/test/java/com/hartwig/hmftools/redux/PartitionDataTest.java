package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.redux.TestUtils.createTestConfig;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.redux.common.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.redux.common.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.UNSET;
import static com.hartwig.hmftools.redux.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createFragment;
import static com.hartwig.hmftools.redux.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.redux.common.CandidateDuplicates;
import com.hartwig.hmftools.redux.common.Fragment;
import com.hartwig.hmftools.redux.common.PartitionData;
import com.hartwig.hmftools.redux.common.PartitionResults;
import com.hartwig.hmftools.redux.common.ResolvedFragmentState;

import org.junit.Test;

public class PartitionDataTest
{
    private final ReadIdGenerator mReadIdGen;
    private final ReduxConfig mConfig;

    private static final String LOCAL_PARTITION_STR = "1_0";

    public PartitionDataTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mConfig = createTestConfig();
    }

    @Test
    public void testBasicFragments()
    {
        PartitionData partitionData = new PartitionData(LOCAL_PARTITION_STR, mConfig);

        List<Fragment> testFragments = createBasicFragments();

        Fragment read = testFragments.get(0);
        Fragment mateRead = testFragments.get(1);
        Fragment supp = testFragments.get(2);

        assertFalse(read.allReadsPresent());

        // test 1: resolved then mate then supp
        List<Fragment> resolvedFragments = Lists.newArrayList(read);
        partitionData.processPrimaryFragments(resolvedFragments, Collections.EMPTY_LIST);
        assertEquals(1, resolvedFragments.size());

        ResolvedFragmentState resolvedState = partitionData.fragmentStatusMap().get(read.id());
        assertNotNull(resolvedState);
        assertEquals(resolvedState.Status, read.status());
        assertFalse(resolvedState.MateReceived);
        assertEquals(1, resolvedState.ExpectedSupplementaries);
        assertEquals(0, resolvedState.ProcessedSupplementaries);

        resolvedFragments = processIncompleteFragment(partitionData, mateRead);
        assertNull(resolvedFragments);
        assertEquals(NONE, mateRead.status());

        resolvedState = partitionData.fragmentStatusMap().get(read.id());
        assertNotNull(resolvedState);
        assertTrue(resolvedState.MateReceived);
        assertEquals(1, resolvedState.ExpectedSupplementaries);
        assertEquals(0, resolvedState.ProcessedSupplementaries);

        resolvedFragments = processIncompleteFragment(partitionData, supp);
        assertNull(resolvedFragments);
        assertEquals(NONE, supp.status());

        assertFalse(partitionData.fragmentStatusMap().containsKey(read.id()));

        partitionData.clearState();

        // test 2: mate then resolved then supp
        testFragments = createBasicFragments();

        read = testFragments.get(0);
        mateRead = testFragments.get(1);
        supp = testFragments.get(2);

        resolvedFragments = processIncompleteFragment(partitionData, mateRead);
        assertNull(resolvedFragments);
        assertEquals(UNSET, mateRead.status());
        Fragment mateFragment = partitionData.incompleteFragmentMap().get(mateRead.id());
        assertTrue(mateFragment != null && mateFragment.reads().contains(mateRead.reads().get(0)));

        resolvedFragments = Lists.newArrayList(read);
        partitionData.processPrimaryFragments(resolvedFragments, Collections.EMPTY_LIST);
        assertEquals(1, resolvedFragments.size());
        assertEquals(2, read.readCount());
        assertFalse(partitionData.incompleteFragmentMap().containsValue(mateRead));

        resolvedState = partitionData.fragmentStatusMap().get(read.id());
        assertNotNull(resolvedState);
        assertEquals(resolvedState.Status, read.status());
        assertTrue(resolvedState.MateReceived);
        assertEquals(1, resolvedState.ExpectedSupplementaries);
        assertEquals(0, resolvedState.ProcessedSupplementaries);

        resolvedFragments = processIncompleteFragment(partitionData, supp);
        assertNull(resolvedFragments);
        assertEquals(NONE, supp.status());

        assertFalse(partitionData.fragmentStatusMap().containsKey(read.id()));

        // test 3: supp then mate then resolved
        testFragments = createBasicFragments();

        read = testFragments.get(0);
        mateRead = testFragments.get(1);
        supp = testFragments.get(2);

        resolvedFragments = processIncompleteFragment(partitionData, supp);
        assertNull(resolvedFragments);
        assertEquals(SUPPLEMENTARY, supp.status());
        Fragment suppFragment = partitionData.incompleteFragmentMap().get(supp.id());
        assertTrue(suppFragment != null && suppFragment.reads().contains(supp.reads().get(0)));
        // assertTrue(partitionData.incompleteFragmentMap().containsValue(supp));

        resolvedFragments = processIncompleteFragment(partitionData, mateRead);
        assertNull(resolvedFragments);
        assertEquals(UNSET, mateRead.status());
        assertEquals(2, suppFragment.readCount());
        assertFalse(partitionData.incompleteFragmentMap().containsValue(mateRead));

        resolvedFragments = Lists.newArrayList(read);
        partitionData.processPrimaryFragments(resolvedFragments, Collections.EMPTY_LIST);
        assertEquals(1, resolvedFragments.size());
        assertEquals(3, read.readCount());
        assertFalse(partitionData.incompleteFragmentMap().containsKey(suppFragment.id()));
        // assertFalse(partitionData.incompleteFragmentMap().containsValue(supp));

        assertFalse(partitionData.fragmentStatusMap().containsKey(read.id()));
    }

    private List<Fragment> createBasicFragments()
    {
        mReadIdGen.reset();

        Fragment read = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read.setStatus(NONE);

        Fragment mateRead = createFragment(read.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null);


        Fragment supp = createFragment(read.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                false, true, new SupplementaryReadData(CHR_1, 2000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        return Lists.newArrayList(read, mateRead, supp);
    }

    @Test
    public void testCandidateFragments()
    {
        PartitionData partitionData = new PartitionData(LOCAL_PARTITION_STR, mConfig);

        // test 1: 2 candidate fragments both waiting on their mates
        List<Fragment> testFragments = createCandidateFragments();
        Fragment read1 = testFragments.get(0);
        Fragment mateRead1 = testFragments.get(1);
        Fragment supp1 = testFragments.get(2);
        Fragment read2 = testFragments.get(3);
        Fragment mateRead2 = testFragments.get(4);

        List<Fragment> resolvedFragments = Lists.newArrayList();
        CandidateDuplicates candidateDuplicates = CandidateDuplicates.from(read1);
        candidateDuplicates.addFragment(read2);
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList(candidateDuplicates);
        partitionData.processPrimaryFragments(resolvedFragments, candidateDuplicatesList);
        assertEquals(0, resolvedFragments.size());

        assertTrue(partitionData.incompleteFragmentMap().containsValue(read1));
        assertTrue(partitionData.incompleteFragmentMap().containsValue(read2));

        // now send through the mate reads in turn
        resolvedFragments = processIncompleteFragment(partitionData, mateRead1);
        assertTrue(resolvedFragments == null || resolvedFragments.isEmpty());
        assertEquals(2, read1.readCount());
        assertTrue(partitionData.incompleteFragmentMap().containsValue(read1));

        resolvedFragments = processIncompleteFragment(partitionData, mateRead2);
        assertNotNull(resolvedFragments);
        assertEquals(2, resolvedFragments.size());
        assertEquals(2, read1.readCount());
        assertFalse(partitionData.incompleteFragmentMap().containsValue(read1));
        assertFalse(partitionData.incompleteFragmentMap().containsValue(read2));

        assertEquals(PRIMARY, read1.status());
        assertEquals(DUPLICATE, read2.status());

        ResolvedFragmentState resolvedState = partitionData.fragmentStatusMap().get(read1.id());
        assertNotNull(resolvedState);
        assertEquals(resolvedState.Status, read1.status());
        assertTrue(resolvedState.MateReceived);
        assertEquals(1, resolvedState.ExpectedSupplementaries);
        assertEquals(0, resolvedState.ProcessedSupplementaries);

        assertFalse(partitionData.fragmentStatusMap().containsKey(read2.id()));

        // finally process the supplementary for the first fragment
        resolvedFragments = processIncompleteFragment(partitionData, supp1);
        assertTrue(resolvedFragments == null || resolvedFragments.isEmpty());
        assertEquals(PRIMARY, supp1.status());
        assertFalse(partitionData.fragmentStatusMap().containsKey(read1.id()));
    }

    private List<Fragment> createCandidateFragments()
    {
        mReadIdGen.reset();

        Fragment read1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setStatus(CANDIDATE);

        Fragment mateRead1 = createFragment(read1.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null);

        Fragment supp1 = createFragment(read1.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                false, true, new SupplementaryReadData(CHR_1, 2000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        Fragment read2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);

        read2.setStatus(CANDIDATE);
        setBaseQualities(read2, DEFAULT_QUAL - 1);

        Fragment mateRead2 = createFragment(read2.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null);

        return Lists.newArrayList(read1, mateRead1, supp1, read2, mateRead2);
    }

    private List<Fragment> processIncompleteFragment(final PartitionData partitionData, final Fragment fragment)
    {
        PartitionResults partitionResults = partitionData.processIncompleteFragment(fragment.reads().get(0));

        if(partitionResults == null)
            return null;

        if(partitionResults.fragmentStatus() != null)
            fragment.setStatus(partitionResults.fragmentStatus());

        return partitionResults.resolvedFragments();
    }

}
