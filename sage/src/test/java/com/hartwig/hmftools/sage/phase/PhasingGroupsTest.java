package com.hartwig.hmftools.sage.phase;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.evidence.PhasedVariantGroup;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.VariantPhaser;

import org.junit.Test;

public class PhasingGroupsTest
{
    private VariantPhaser mPhaser;
    private int mNextReadCounterId;

    public PhasingGroupsTest()
    {
        mPhaser = new VariantPhaser(new PhaseSetCounter());
        mPhaser.initialise(new ChrBaseRegion("1", 0, 200), true, false);

        mNextReadCounterId = 0;
    }

    @Test
    public void testGroupCreations()
    {
        // test that groups are created distinctly where not exact matches, ordered by start position
        ReadContextCounter rc1 = createReadCounter(10);
        ReadContextCounter rc2 = createReadCounter(10);
        ReadContextCounter rc3 = createReadCounter(20);
        ReadContextCounter rc4 = createReadCounter(20);
        ReadContextCounter rc5 = createReadCounter(30);
        ReadContextCounter rc6 = createReadCounter(40);

        List<ReadContextCounter> posCounters = Lists.newArrayList(rc1);
        List<ReadContextCounter> negCounters = Lists.newArrayList(rc2);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(1, mPhaser.getPhasedGroups().size());

        // negatives can diff
        posCounters = Lists.newArrayList(rc1);
        negCounters = Lists.newArrayList(rc2, rc3);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(1, mPhaser.getPhasedGroups().size());

        PhasedVariantGroup group = findGroup(posCounters);
        assertNotNull(group);
        assertEquals(2, group.ReadCount);

        // positive cannot differ
        posCounters = Lists.newArrayList(rc1, rc3);
        negCounters = Lists.newArrayList(rc2);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(2, mPhaser.getPhasedGroups().size());
        assertTrue(checkGroupsAreOrdered());

        // groups are added in order
        posCounters = Lists.newArrayList(rc5, rc6);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(3, mPhaser.getPhasedGroups().size());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc3, rc5);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(4, mPhaser.getPhasedGroups().size());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc2, rc4);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(5, mPhaser.getPhasedGroups().size());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc6);
        negCounters = Lists.newArrayList(rc3);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(6, mPhaser.getPhasedGroups().size());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc4, rc5);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(7, mPhaser.getPhasedGroups().size());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc2, rc5);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(8, mPhaser.getPhasedGroups().size());
        assertTrue(checkGroupsAreOrdered());

        // test matches on positives only and merges of negatives
        posCounters = Lists.newArrayList(rc6);
        negCounters = Lists.newArrayList(rc3, rc4, rc5);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(8, mPhaser.getPhasedGroups().size());
        group = findGroup(posCounters);
        assertNotNull(group);
        assertEquals(2, group.ReadCount);
        assertEquals(3, group.NegativeReadCounters.size());

        posCounters = Lists.newArrayList(rc4, rc5);
        negCounters = Lists.newArrayList(rc1, rc6);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(8, mPhaser.getPhasedGroups().size());
        group = findGroup(posCounters);
        assertNotNull(group);
        assertEquals(2, group.ReadCount);
        assertEquals(2, group.NegativeReadCounters.size());
    }

    @Test
    public void testGroupMerging()
    {
        // scenarios:
        // subsets allocated pro-rate to super group(s), test with subsets of subsets of super groups
        // single-variant groups dropped if not present in another group
        // merge of groups with common subset
        // cannot merge if have conflicting negatives

        ReadContextCounter rc0 = createReadCounter(10);
        ReadContextCounter rc1 = createReadCounter(10);
        ReadContextCounter rc2 = createReadCounter(20);
        ReadContextCounter rc3 = createReadCounter(20);
        ReadContextCounter rc4 = createReadCounter(30);
        ReadContextCounter rc5 = createReadCounter(40);
        ReadContextCounter rc6 = createReadCounter(50);
        ReadContextCounter rc7 = createReadCounter(60);

        List<ReadContextCounter> posCounters = Lists.newArrayList(rc0);
        List<ReadContextCounter> negCounters = Lists.newArrayList(rc1);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(1, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group0 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        // not super group of 1 since has conflicting negative
        posCounters = Lists.newArrayList(rc0, rc1);
        negCounters = Lists.newArrayList(rc2);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(2, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group1 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        // super group of 1
        posCounters = Lists.newArrayList(rc0, rc2);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(3, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group2 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        // another super group of 1, and subgroup of 4, 5 and 6
        posCounters = Lists.newArrayList(rc0, rc3);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(4, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group3 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        // super group of 4
        posCounters = Lists.newArrayList(rc0, rc3, rc4);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(5, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group4 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        // another super group of 4, will merge with group 5 since has common subset and no conflicts
        posCounters = Lists.newArrayList(rc0, rc3, rc5);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(6, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group5 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        // another super group of 4 but cannot be merged with common subsets
        posCounters = Lists.newArrayList(rc0, rc3, rc7);
        negCounters = Lists.newArrayList(rc4, rc5);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(7, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group6 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        // single group not present in others will be dropped
        posCounters = Lists.newArrayList(rc6);
        negCounters = Lists.newArrayList(rc4, rc5);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(8, mPhaser.getPhasedGroups().size());
        PhasedVariantGroup group7 = mPhaser.getPhasedGroups().get(mPhaser.getPhasedGroups().size() - 1);

        mPhaser.mergeGroups();

        /* when -ves are not copied from subsets to supersets
        assertEquals(2, mPhaser.getPhasedGroups().size());

        PhasedVariantGroup group = findGroup(Lists.newArrayList(rc0, rc1, rc3, rc4, rc5));
        assertNotNull(group);
        assertEquals(3, group.ReadCount);
        assertEquals(1.2, group.AllocatedReadCount, 0.1);

        group = findGroup(Lists.newArrayList(rc0, rc2, rc3, rc7));
        assertNotNull(group);
        assertEquals(2, group.ReadCount);
        assertEquals(0.8, group.AllocatedReadCount, 0.1);
        */

        assertEquals(3, mPhaser.getPhasedGroups().size());

        PhasedVariantGroup group = findGroup(Lists.newArrayList(rc0, rc3, rc4, rc5));
        assertNotNull(group);
        assertEquals(2, group.ReadCount);
        assertEquals(1.2, group.AllocatedReadCount, 0.1);

        group = findGroup(Lists.newArrayList(rc0, rc1));
        assertNotNull(group);
        assertEquals(1, group.ReadCount);
        assertEquals(0, group.AllocatedReadCount, 0.1);

        group = findGroup(Lists.newArrayList(rc0, rc2, rc3, rc7));
        assertNotNull(group);
        assertEquals(2, group.ReadCount);
        assertEquals(0.8, group.AllocatedReadCount, 0.1);
    }

    private PhasedVariantGroup findGroup(final List<ReadContextCounter> posCounters)
    {
        int minPos = posCounters.get(0).position();
        int maxPos = posCounters.get(posCounters.size() - 1).position();
        return mPhaser.getPhasedGroups().stream().filter(x -> x.matches(minPos, maxPos, posCounters)).findFirst().orElse(null);
    }

    private boolean checkGroupsAreOrdered()
    {
        for(int i = 0; i < mPhaser.getPhasedGroups().size() - 1; ++i)
        {
            PhasedVariantGroup current = mPhaser.getPhasedGroups().get(i);
            PhasedVariantGroup next = mPhaser.getPhasedGroups().get(i + 1);

            if(next.minVariantPos() < current.minVariantPos())
                return false;
        }

        return true;
    }

    private ReadContextCounter createReadCounter(int position)
    {
        return createReadCounter(mNextReadCounterId++, position);
    }

    private ReadContextCounter createReadCounter(int id, int position)
    {
        VariantHotspot variant = ImmutableVariantHotspotImpl.builder()
                .chromosome("1")
                .position(position)
                .ref("A")
                .alt("C").build();

        IndexedBases indexBases = new IndexedBases(position, 10, "ACGTACGTACGT".getBytes());
        ReadContext readContext = new ReadContext(position, "", 0, "", indexBases, false);

        return new ReadContextCounter(id, variant, readContext, VariantTier.LOW_CONFIDENCE,
                100, 1, false);
    }
}
