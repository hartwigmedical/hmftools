package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragment;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.junit.Test;

public class LocalGroupCombinerTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MarkDupsConfig mConfig;
    private final RecordWriter mWriter;

    private static final String LOCAL_PARTITION_STR = "1_0";
    private static final BaseRegion LOCAL_PARTITION = new BaseRegion(1, 1000000);

    public LocalGroupCombinerTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mConfig = new MarkDupsConfig();
        mWriter = new RecordWriter(mConfig);
    }

    @Test
    public void testLocalFragmentBasic()
    {
        GroupCombiner localGroupCombiner = new GroupCombiner(mWriter, true, false);

        List<Fragment> testFragments = createBasicFragments();

        Fragment read1 = testFragments.get(0);
        read1.setRemotePartitions(LOCAL_PARTITION);

        List<Fragment> resolvedFragments = Lists.newArrayList(read1);
        List<CandidateDuplicates> incompletePositionFragments = Lists.newArrayList();
        List<Fragment> supplementaries = Lists.newArrayList();

        // test 1: resolved then supps then unclear
        supplementaries.add(testFragments.get(2));
        supplementaries.add(testFragments.get(3));
        incompletePositionFragments.add(new CandidateDuplicates(testFragments.get(1)));

        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, resolvedFragments, Collections.emptyList(), Collections.emptyList());
        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), Collections.emptyList(), supplementaries);
        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), incompletePositionFragments, Collections.emptyList());

        List<Fragment> targets = Lists.newArrayList(testFragments.get(1), testFragments.get(2), testFragments.get(3));
        targets.forEach(x -> assertTrue(x.readsWritten()));
        targets.forEach(x -> assertEquals(NONE, x.status()));

        localGroupCombiner.reset();

        // reads have been merged and status changed, need to recreate
        testFragments = createBasicFragments();

        supplementaries.clear();
        supplementaries.add(testFragments.get(2));
        supplementaries.add(testFragments.get(3));
        Fragment read2 = testFragments.get(1);
        incompletePositionFragments.add(new CandidateDuplicates(read2));
        targets = Lists.newArrayList(testFragments.get(1), testFragments.get(2), testFragments.get(3));

        // test 2: supps then unclear then resolved
        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), Collections.emptyList(), supplementaries);
        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), incompletePositionFragments, Collections.emptyList());
        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, resolvedFragments, Collections.emptyList(), Collections.emptyList());

        assertTrue(read2.readsWritten());
        assertEquals(NONE, read2.status());
        assertEquals(3, read2.readCount()); // supplementaries were merged

        localGroupCombiner.reset();

        // test 3: unclear then supp then resolved
        testFragments = createBasicFragments();

        supplementaries.clear();
        supplementaries.add(testFragments.get(2));
        supplementaries.add(testFragments.get(3));
        read2 = testFragments.get(1);
        incompletePositionFragments.add(new CandidateDuplicates(read2));
        targets = Lists.newArrayList(testFragments.get(1), testFragments.get(2), testFragments.get(3));

        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), incompletePositionFragments, Collections.emptyList());
        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), Collections.emptyList(), supplementaries);
        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, resolvedFragments, Collections.emptyList(), Collections.emptyList());

        assertTrue(read2.readsWritten());
        assertEquals(NONE, read2.status());
        assertEquals(3, read2.readCount()); // supplementaries were merged
    }

    private List<Fragment> createBasicFragments()
    {
        mReadIdGen.reset();

        Fragment read1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setStatus(NONE);

        Fragment read2 = createFragment(read1.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setStatus(UNCLEAR);

        Fragment supp1 = createFragment(read1.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                false, true, new SupplementaryReadData(CHR_1, 2000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        Fragment supp2 = createFragment(read1.id(), CHR_1, 1000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, true, new SupplementaryReadData(CHR_1, 2000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        return Lists.newArrayList(read1, read2, supp1, supp2);
    }

    @Test
    public void testLocalFragmentSplitReads()
    {
        GroupCombiner localGroupCombiner = new GroupCombiner(mWriter, true, false);

        // a collection of read pairs, all unclear and split across calls
        mReadIdGen.reset();

        // 3 initial reads at the same position, unclear candidate duplicates
        Fragment read1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, null);

        // has a supp that will be passed in later on
        Fragment read2 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, false, new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));
        setBaseQualities(read2, DEFAULT_QUAL - 1);

        Fragment suppRead2 = createFragment(read2.id(), CHR_1, 1000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 201,
                false, true, new SupplementaryReadData(CHR_1, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        Fragment read3 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 202,
                false, false, null);
        setBaseQualities(read3, DEFAULT_QUAL - 1);

        List<Fragment> positionFragmentsList = Lists.newArrayList(read1, read2, read3);
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> incompletePositionFragments = Lists.newArrayList();

        classifyFragments(positionFragmentsList, resolvedFragments, incompletePositionFragments);
        assertEquals(1, incompletePositionFragments.size());
        assertEquals(3, incompletePositionFragments.get(0).Fragments.size());

        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), incompletePositionFragments, Collections.emptyList());
        assertEquals(1, localGroupCombiner.getPartitionCache(LOCAL_PARTITION_STR).CandidateDuplicatesMap.size());

        // now process their mates one by one

        // first 2 are also UNCLEAR and will resolve as duplicate with their mates
        Fragment readMate1 = createFragment(read1.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                true, false, null);

        Fragment readMate2 = createFragment(read2.id(), CHR_1, 201, TEST_READ_BASES.substring(0, 99), "99M", CHR_1, 100,
                true, false, null);
        setBaseQualities(readMate2, DEFAULT_QUAL - 1);

        positionFragmentsList = Lists.newArrayList(readMate1, readMate2);
        resolvedFragments.clear();
        incompletePositionFragments.clear();

        classifyFragments(positionFragmentsList, resolvedFragments, incompletePositionFragments);
        assertEquals(1, incompletePositionFragments.size());
        assertEquals(2, incompletePositionFragments.get(0).Fragments.size());

        localGroupCombiner.processPartitionFragments(LOCAL_PARTITION_STR, Collections.emptyList(), incompletePositionFragments, Collections.emptyList());

        assertTrue(readMate1.readsWritten());
        assertEquals(PRIMARY, readMate1.status());
        assertTrue(readMate2.readsWritten());
        assertEquals(DUPLICATE, readMate2.status());

        // 3rd fragment now able to be classified as NONE since no other candidate dups exist at this position
        assertTrue(read3.readsWritten());
        assertEquals(NONE, read3.status());

        assertTrue(localGroupCombiner.getPartitionCache(LOCAL_PARTITION_STR).IncompleteFragments.isEmpty());

        // finally pass in the supplementary which should pick up the resolved status
        localGroupCombiner.localSupplementary(suppRead2, LOCAL_PARTITION_STR);

        assertTrue(suppRead2.readsWritten());
        assertEquals(DUPLICATE, suppRead2.status());
    }
}
