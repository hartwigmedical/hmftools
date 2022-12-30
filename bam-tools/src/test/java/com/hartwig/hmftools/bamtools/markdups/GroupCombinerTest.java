package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.bamtools.markdups.TestUtils.createFragment;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.junit.Test;

public class GroupCombinerTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MarkDupsConfig mConfig;
    private final RecordWriter mWriter;

    private static final String LOCAL_PARTITION_STR = "1_0";
    private static final BaseRegion LOCAL_PARTITION = new BaseRegion(1, 1000000);

    public GroupCombinerTest()
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
        Fragment supp1 = testFragments.get(0);
        read1.setRemotePartitions(LOCAL_PARTITION);
        assertTrue(read1.hasRemotePartitions());

        // test 1: resolved then supps
        localGroupCombiner.localResolvedFragment(LOCAL_PARTITION_STR, read1);
        localGroupCombiner.localSupplementary(LOCAL_PARTITION_STR, supp1);

        assertTrue(supp1.readsWritten());
        assertEquals(NONE, supp1.status());

        localGroupCombiner.reset();

        // reads have been merged and status changed, need to recreate
        testFragments = createBasicFragments();

        read1 = testFragments.get(0);
        supp1 = testFragments.get(0);

        // test 2: supps then resolved
        localGroupCombiner.localSupplementary(LOCAL_PARTITION_STR, supp1);

        read1.setRemotePartitions(LOCAL_PARTITION);
        localGroupCombiner.localResolvedFragment(LOCAL_PARTITION_STR, read1);

        assertTrue(supp1.readsWritten());
        assertEquals(NONE, supp1.status());
    }

    private List<Fragment> createBasicFragments()
    {
        mReadIdGen.reset();

        Fragment read1 = createFragment(mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 200,
                false, false, new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setStatus(NONE);

        Fragment supp1 = createFragment(read1.id(), CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                false, true, new SupplementaryReadData(CHR_1, 2000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        return Lists.newArrayList(read1, supp1);
    }

}
