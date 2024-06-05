package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.redux.TestUtils.createFragment;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.redux.umi.UmiConfig.extractUmiIdFromReadId;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.umi.UmiGroupBuilder.buildUmiGroups;
import static com.hartwig.hmftools.redux.umi.UmiGroupBuilder.hasDuplexUmiMatch;
import static com.hartwig.hmftools.redux.umi.UmiUtils.exceedsUmiIdDiff;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.redux.common.Fragment;
import com.hartwig.hmftools.redux.umi.UmiConfig;
import com.hartwig.hmftools.redux.common.DuplicateGroup;

import org.junit.Test;

public class UmiGroupsTest
{
    public static final String FIXED_READ_ID = "123:ABC:1:4455:";

    private static final UmiConfig UMI_CONFIG = new UmiConfig(true, false, "", false);

    @Test
    public void testUmiUtils()
    {
        String umiId = "TATCGC";
        String readId = FIXED_READ_ID + umiId;
        assertEquals(umiId, extractUmiIdFromReadId(readId));

        short umiLength = (short)umiId.length();
        assertEquals(umiId, UMI_CONFIG.extractUmiId(readId, umiLength));

        String readId1 = FIXED_READ_ID + "TATCGC";
        String readId2 = FIXED_READ_ID + "TATCGG";
        String readId3 = FIXED_READ_ID + "TATGGG";
        String readId4 = FIXED_READ_ID + "AAACGG";
        String readId5 = FIXED_READ_ID + "AAAGGG";
        String readId6 = FIXED_READ_ID + "AAAGGGT";

        assertFalse(exceedsUmiIdDiff(readId1, readId2));
        assertFalse(exceedsUmiIdDiff(readId2, readId3));
        assertTrue(exceedsUmiIdDiff(readId1, readId3));
        assertTrue(exceedsUmiIdDiff(readId2, readId4));
        assertFalse(exceedsUmiIdDiff(readId4, readId5));
        assertTrue(exceedsUmiIdDiff(readId5, readId6));
    }

    @Test
    public void testDuplexUmis()
    {
        UmiConfig umiConfig = new UmiConfig(true, true, String.valueOf(DEFAULT_DUPLEX_UMI_DELIM), false);

        String umiId1 = "TATCGC_AAGTCG";
        assertFalse(hasDuplexUmiMatch(umiId1, umiId1, umiConfig.DuplexDelim, umiConfig.PermittedBaseDiff));

        String umiId2 = "AAGTCC_TATCGG";

        assertTrue(hasDuplexUmiMatch(umiId1, umiId2, umiConfig.DuplexDelim, umiConfig.PermittedBaseDiff));

        // too manu different base
        umiId2 = "AAGCCC_TATCGC";

        assertFalse(hasDuplexUmiMatch(umiId1, umiId2, umiConfig.DuplexDelim, umiConfig.PermittedBaseDiff));
    }

    @Test
    public void testUmiGroupAssignment()
    {
        Fragment frag1 = createFragment(FIXED_READ_ID + "TATCGC", CHR_1, 100);
        Fragment frag2 = createFragment(FIXED_READ_ID + "TATCGC", CHR_1, 100);
        Fragment frag3 = createFragment(FIXED_READ_ID + "TATCCC", CHR_1, 100);
        Fragment frag4 = createFragment(FIXED_READ_ID + "AATGGG", CHR_1, 100);
        Fragment frag5 = createFragment(FIXED_READ_ID + "TATCGC", CHR_1, 100);
        Fragment frag6 = createFragment(FIXED_READ_ID + "AATGGG", CHR_1, 100);

        List<Fragment> fragments = Lists.newArrayList(frag1, frag2, frag3, frag4, frag5, frag6);
        List<DuplicateGroup> groups = buildUmiGroups(fragments, UMI_CONFIG);
        assertEquals(2, groups.size());

        DuplicateGroup group = groups.stream().filter(x -> x.fragments().contains(frag1)).findFirst().orElse(null);
        assertTrue(group.fragments().contains(frag2));
        assertTrue(group.fragments().contains(frag3));
        assertTrue(group.fragments().contains(frag5));

        group = groups.stream().filter(x -> x.fragments().contains(frag4)).findFirst().orElse(null);
        assertTrue(group.fragments().contains(frag6));

        // check that larger groups aren't merged
        Fragment frag11 = createFragment(FIXED_READ_ID + "TTTCGT", CHR_1, 100);
        Fragment frag12 = createFragment(FIXED_READ_ID + "TTCCGT", CHR_1, 100);
        Fragment frag13 = createFragment(FIXED_READ_ID + "TTACGT", CHR_1, 100);
        Fragment frag14 = createFragment(FIXED_READ_ID + "TTACAT", CHR_1, 100);
        Fragment frag15 = createFragment(FIXED_READ_ID + "TAAAAT", CHR_1, 100);
        Fragment frag16 = createFragment(FIXED_READ_ID + "TTACAG", CHR_1, 100);

        // duplicate some to build bigger UMI groups
        fragments = Lists.newArrayList(
                frag11, frag11,
                frag12, frag12, frag12,
                frag13, frag13, frag13,frag13, frag13, frag13, frag13, frag13,
                frag14, frag14, frag14, frag14,
                frag15, frag15, frag15, frag15, frag15,
                frag16);

        groups = buildUmiGroups(fragments, UMI_CONFIG);
        assertEquals(2, groups.size());
        assertEquals(fragments.size(), groups.stream().mapToInt(x -> x.fragmentCount()).sum());

        group = groups.stream().filter(x -> x.fragments().contains(frag11)).findFirst().orElse(null);
        assertEquals(18, group.fragmentCount());
        assertTrue(group.fragments().contains(frag12));
        assertTrue(group.fragments().contains(frag13));
        assertTrue(group.fragments().contains(frag14));
        assertTrue(group.fragments().contains(frag16));

        group = groups.stream().filter(x -> x.fragments().contains(frag15)).findFirst().orElse(null);
        assertEquals(5, group.fragmentCount());
    }

    @Test
    public void testUmiGroupAssignment2()
    {
        // test a final merge of groups with 2 bases difference
        Fragment frag1 = createFragment(FIXED_READ_ID + "TTAAGG", CHR_1, 100);
        Fragment frag2 = createFragment(FIXED_READ_ID + "TTAAGC", CHR_1, 100);
        Fragment frag3 = createFragment(FIXED_READ_ID + "TTAATT", CHR_1, 100);
        Fragment frag4 = createFragment(FIXED_READ_ID + "GGGATT", CHR_1, 100);
        Fragment frag5 = createFragment(FIXED_READ_ID + "GGGAGG", CHR_1, 100);
        Fragment frag6 = createFragment(FIXED_READ_ID + "CCGAGC", CHR_1, 100);

        List<Fragment> fragments = Lists.newArrayList(frag1, frag2, frag3, frag4, frag5, frag6);
        List<DuplicateGroup> groups = buildUmiGroups(fragments, UMI_CONFIG);
        assertEquals(3, groups.size());
    }

    @Test
    public void testUmiGroupAssignment3()
    {
        // test a final merge of groups with 4 bases difference if the groups are large enough
        String umi1 = "TTAAGG";
        String umi2 = "TTAGGA"; // 2-base diff from #1
        String umi3 = "CCGGTT"; // unrelated
        String umi4 = "TTCCCC"; // 4-base diff from #1 and #2

        List<Fragment> fragments = Lists.newArrayList();

        fragments.add(createFragment(FIXED_READ_ID + umi1, CHR_1, 100));
        fragments.add(createFragment(FIXED_READ_ID + umi2, CHR_1, 100));
        fragments.add(createFragment(FIXED_READ_ID + umi3, CHR_1, 100));

        for(int i = 0; i < 51; ++i)
        {
            fragments.add(createFragment(FIXED_READ_ID + umi4, CHR_1, 100));
        }

        List<DuplicateGroup> groups = buildUmiGroups(fragments, UMI_CONFIG);
        assertEquals(2, groups.size());
    }

    @Test
    public void testDefinedUmis()
    {
        UmiConfig umiConfig = new UmiConfig(true, false, "", false);
        String definedUmi1 = "AAAGGG";
        String definedUmi2 = "TTTAAA";
        String definedUmi3 = "CCCAAA";
        umiConfig.addDefinedUmis(Sets.newHashSet(definedUmi1, definedUmi2, definedUmi3));

        Fragment frag1 = createFragment(FIXED_READ_ID + definedUmi1, CHR_1, 100);
        Fragment frag2 = createFragment(FIXED_READ_ID + "AAAGGC", CHR_1, 100);
        Fragment frag3 = createFragment(FIXED_READ_ID + definedUmi1, CHR_1, 100);
        Fragment frag4 = createFragment(FIXED_READ_ID + definedUmi2, CHR_1, 100);
        Fragment frag5 = createFragment(FIXED_READ_ID + "TTTAAC", CHR_1, 100);
        Fragment frag6 = createFragment(FIXED_READ_ID + definedUmi3, CHR_1, 100);

        List<Fragment> fragments = Lists.newArrayList(frag1, frag2, frag3, frag4, frag5, frag6);
        List<DuplicateGroup> groups = buildUmiGroups(fragments, umiConfig);
        assertEquals(3, groups.size());
    }

    /*
    @Test
    public void testPerfUmiIdExtraction()
    {
        String umiId = "TATCGC";
        String readId = FIXED_READ_ID + umiId;
        short umiLength = (short) umiId.length();

        int testIterations = 1000000;
        long startTimeMs = System.currentTimeMillis();

        for(int i = 0; i < testIterations; ++i)
        {
            extractUmiId(readId);
        }

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;

        BM_LOGGER.info("slow time: {}", timeTakenMs);

        startTimeMs = System.currentTimeMillis();

        for(int i = 0; i < testIterations; ++i)
        {
            extractUmiId(readId, umiLength);
        }

        long timeTakenMs2 = System.currentTimeMillis() - startTimeMs;

        BM_LOGGER.info("fast time: {}", timeTakenMs2);
    }
    */
}
