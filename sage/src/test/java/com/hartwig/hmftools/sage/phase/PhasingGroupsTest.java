package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.config.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.phase.VariantPhaser.mergeByExtension;
import static com.hartwig.hmftools.sage.phase.VariantPhaser.mergeMatching;
import static com.hartwig.hmftools.sage.phase.VariantPhaser.mergeUninformative;
import static com.hartwig.hmftools.sage.phase.VariantPhaser.removeUninformativeLps;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.junit.Test;

public class PhasingGroupsTest
{
    private VariantPhaser mPhaser;
    private int mNextReadCounterId;

    public PhasingGroupsTest()
    {
        mPhaser = new VariantPhaser(new PhaseSetCounter());
        mPhaser.initialise(new ChrBaseRegion("1", 0, 200), false);

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
        assertEquals(1, getPhasedGroupCount());

        // negatives cannot diff
        posCounters = Lists.newArrayList(rc1);
        negCounters = Lists.newArrayList(rc2, rc3);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(2, getPhasedGroupCount());

        PhasedVariantGroup group = findGroup(posCounters);
        assertNotNull(group);
        assertEquals(1, group.ReadCount);

        // positive cannot differ
        posCounters = Lists.newArrayList(rc1, rc3);
        negCounters = Lists.newArrayList(rc2);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(3, getPhasedGroupCount());
        assertTrue(checkGroupsAreOrdered());

        // groups are added in order
        posCounters = Lists.newArrayList(rc5, rc6);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(4, getPhasedGroupCount());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc3, rc5);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(5, getPhasedGroupCount());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc2, rc4);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(6, getPhasedGroupCount());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc6);
        negCounters = Lists.newArrayList(rc3);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(7, getPhasedGroupCount());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc4, rc5);
        negCounters = Lists.newArrayList();

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(8, getPhasedGroupCount());
        assertTrue(checkGroupsAreOrdered());

        posCounters = Lists.newArrayList(rc2, rc5);

        mPhaser.registeredPhasedVariants(posCounters, negCounters);
        assertEquals(9, getPhasedGroupCount());
        assertTrue(checkGroupsAreOrdered());
    }

    private int getPhasedGroupCount()
    {
        return mPhaser.getPhasedCollections().stream().mapToInt(x -> x.groups().size()).sum();
    }

    @Test
    public void testMergeMatching()
    {
        // scenarios:
        // remove groups which only have their +ves in other groups with the same +ves
        // single-variant groups dropped if not present in another group

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

        // matches group 0, will be merged
        posCounters = Lists.newArrayList(rc0);
        negCounters = Lists.newArrayList(rc1, rc2);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // also matches group 0, will be merged
        posCounters = Lists.newArrayList(rc0);
        negCounters = Lists.newArrayList(rc2);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // distinct group and will be a subset of another
        posCounters = Lists.newArrayList(rc1, rc2);
        negCounters = Lists.newArrayList(rc3);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // superset of 3
        posCounters = Lists.newArrayList(rc1, rc2, rc4);
        negCounters = Lists.newArrayList(rc5);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // matches group 4, will be merged
        posCounters = Lists.newArrayList(rc1, rc2, rc4);
        negCounters = Lists.newArrayList(rc6);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // finally 2 groups not part of any others
        posCounters = Lists.newArrayList(rc3);
        negCounters = Lists.newArrayList(rc6);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        posCounters = Lists.newArrayList(rc6);
        negCounters = Lists.newArrayList(rc7);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        List<PhasedVariantGroup> groups = Lists.newArrayList(getPhasedGroups());

        assertEquals(8, groups.size());

        mergeMatching(groups);
        assertEquals(3, groups.size());

        PhasedVariantGroup group = findGroup(Lists.newArrayList(rc0));
        assertNotNull(group);
        assertEquals(3, group.ReadCount);
        assertEquals(0, group.AllocatedReadCount, 0.1);

        group = findGroup(Lists.newArrayList(rc1, rc2));
        assertNotNull(group);
        assertEquals(1, group.ReadCount);

        group = findGroup(Lists.newArrayList(rc1, rc2, rc4));
        assertNotNull(group);
        assertEquals(2, group.ReadCount);
    }

    @Test
    public void testMergByExtension()
    {
        // merge any group with a common subset of +ves and -ves if it can only be extended in one direction
        ReadContextCounter rc0 = createReadCounter(10);
        ReadContextCounter rc1 = createReadCounter(20);
        ReadContextCounter rc2 = createReadCounter(30);
        ReadContextCounter rc3 = createReadCounter(40);
        ReadContextCounter rc4 = createReadCounter(50);
        ReadContextCounter rc5 = createReadCounter(60);
        ReadContextCounter rc6 = createReadCounter(70);
        // ReadContextCounter rc7 = createReadCounter(80);

        // 2 sub-groups and 1 extending in either direction, results in a single group
        List<ReadContextCounter> posCounters = Lists.newArrayList(rc2, rc5);
        List<ReadContextCounter> negCounters = Lists.newArrayList(rc3, rc4);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // another sub-group
        posCounters = Lists.newArrayList(rc2, rc5);
        negCounters = Lists.newArrayList(rc4);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // extends down
        posCounters = Lists.newArrayList(rc0, rc1, rc2, rc5);
        negCounters = Lists.newArrayList();
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // extends up
        posCounters = Lists.newArrayList(rc2, rc5, rc6);
        negCounters = Lists.newArrayList();
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        List<PhasedVariantGroup> groups = Lists.newArrayList(getPhasedGroups());

        assertEquals(4, groups.size());

        mergeByExtension(groups);
        assertEquals(2, groups.size());

        PhasedVariantGroup group = findGroup(Lists.newArrayList(rc0, rc1, rc2, rc5, rc6));
        assertNotNull(group);
        assertEquals(1, group.ReadCount);
        assertEquals(1, group.AllocatedReadCount, 0.1);

        mergeUninformative(groups);
        assertEquals(1, groups.size());


        // test again but with multiple options up and down preventing a merge
        groups.clear();
        mPhaser.initialise(mPhaser.region(), false);

        // a sub-group
        posCounters = Lists.newArrayList(rc2, rc3, rc4);
        negCounters = Lists.newArrayList(rc4);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // conflicting sub-group
        posCounters = Lists.newArrayList(rc2, rc4);
        negCounters = Lists.newArrayList(rc3);
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // extends down
        posCounters = Lists.newArrayList(rc1, rc2, rc3, rc4);
        negCounters = Lists.newArrayList();
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // extends down again differently
        posCounters = Lists.newArrayList(rc0, rc2, rc3, rc4);
        negCounters = Lists.newArrayList();
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // extends up
        posCounters = Lists.newArrayList(rc2, rc3, rc4, rc5);
        negCounters = Lists.newArrayList();
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        // extends up differently
        posCounters = Lists.newArrayList(rc2, rc3, rc4, rc6);
        negCounters = Lists.newArrayList();
        mPhaser.registeredPhasedVariants(posCounters, negCounters);

        groups = Lists.newArrayList(getPhasedGroups());

        assertEquals(6, groups.size());

        mergeByExtension(groups);

        assertEquals(6, groups.size());
    }

    @Test
    public void testRemoveUninformativeLps()
    {
        Set<Integer> passingPhaseSets = Sets.newHashSet();
        List<SageVariant> variants = Lists.newArrayList();

        // LPS 1 and 2 both have passing var 1
        // LPS 2 has same passing variant as LPS 1 and is uninformative

        // LPS 3, 4 and 5 all have passing vars 4 & 5
        // LPS 4 & 5 has same passing variants as LPS 3 and are uninformative

        // LPS 6 has passing variants 1 and 4 but no overlap so remains

        SageVariant var1 = createVariant(11);
        var1.tumorReadCounters().get(0).addLocalPhaseSet(1, 10, 0);
        var1.tumorReadCounters().get(0).addLocalPhaseSet(2, 1, 0);
        var1.tumorReadCounters().get(0).addLocalPhaseSet(6, 1, 0);
        variants.add(var1);

        SageVariant var2 = createVariant(12);
        var2.tumorReadCounters().get(0).addLocalPhaseSet(1, 10, 0);
        var2.filters().add(MAX_GERMLINE_VAF.toString());
        variants.add(var2);

        // a second group, uninformative vs both 1 and 3
        SageVariant var3 = createVariant(13);
        var3.tumorReadCounters().get(0).addLocalPhaseSet(2, 1, 0);
        var3.filters().add(MAX_GERMLINE_VAF.toString());
        variants.add(var3);

        SageVariant var4 = createVariant(14);
        var4.tumorReadCounters().get(0).addLocalPhaseSet(3, 10, 0);
        var4.tumorReadCounters().get(0).addLocalPhaseSet(4, 5, 0);
        var4.tumorReadCounters().get(0).addLocalPhaseSet(5, 1, 0);
        var4.tumorReadCounters().get(0).addLocalPhaseSet(6, 1, 0);
        variants.add(var4);

        SageVariant var5 = createVariant(15);
        var5.tumorReadCounters().get(0).addLocalPhaseSet(3, 10, 0);
        var5.tumorReadCounters().get(0).addLocalPhaseSet(4, 5, 0);
        var5.tumorReadCounters().get(0).addLocalPhaseSet(5, 1, 0);

        // a group which is a subset and with less < 25% RCs
        var5.tumorReadCounters().get(0).addLocalPhaseSet(7, 1, 0);

        variants.add(var5);

        SageVariant var6 = createVariant(16);
        var6.tumorReadCounters().get(0).addLocalPhaseSet(3, 10, 0);
        var6.filters().add(MAX_GERMLINE_VAF.toString());
        variants.add(var6);

        SageVariant var7 = createVariant(17);
        var7.tumorReadCounters().get(0).addLocalPhaseSet(4, 5, 0);
        var7.filters().add(MAX_GERMLINE_VAF.toString());
        variants.add(var7);

        SageVariant var8 = createVariant(18);
        var8.tumorReadCounters().get(0).addLocalPhaseSet(5, 1, 0);
        var8.tumorReadCounters().get(0).addLocalPhaseSet(6, 10, 0);
        var8.filters().add(MAX_GERMLINE_VAF.toString());
        variants.add(var8);

        SageVariant var9 = createVariant(19);
        var9.tumorReadCounters().get(0).addLocalPhaseSet(7, 1, 0);
        var9.filters().add(MAX_GERMLINE_VAF.toString());
        variants.add(var9);

        variants.stream().filter(x -> x.isPassing()).forEach(x -> x.localPhaseSets().stream().forEach(y -> passingPhaseSets.add(y)));
        assertEquals(7, passingPhaseSets.size());

        removeUninformativeLps(variants, passingPhaseSets);

        assertTrue(var1.hasMatchingLps(1));

        assertFalse(var1.hasMatchingLps(2));
        assertFalse(var3.hasMatchingLps(2));

        assertTrue(var4.hasMatchingLps(3));
        assertTrue(var5.hasMatchingLps(3));
        assertTrue(var6.hasMatchingLps(3));

        assertFalse(var4.hasMatchingLps(4));
        assertFalse(var5.hasMatchingLps(4));
        assertFalse(var7.hasMatchingLps(4));

        assertFalse(var4.hasMatchingLps(5));
        assertFalse(var5.hasMatchingLps(5));
        assertFalse(var9.hasMatchingLps(5));

        assertTrue(var1.hasMatchingLps(6));
        assertTrue(var4.hasMatchingLps(6));

        assertFalse(var5.hasMatchingLps(7));
        assertFalse(var9.hasLocalPhaseSets());
    }

    private List<PhasedVariantGroup> getPhasedGroups()
    {
        List<PhasedVariantGroup> phasedGroups = Lists.newArrayList();
        mPhaser.getPhasedCollections().forEach(x -> phasedGroups.addAll(x.groups()));
        return phasedGroups;
    }

    private PhasedVariantGroup findGroup(final List<ReadContextCounter> posCounters)
    {
        return getPhasedGroups().stream()
                .filter(x -> x.PositiveReadCounters.size() == posCounters.size())
                .filter(x -> x.PositiveReadCounters.stream().allMatch(y -> posCounters.contains(y)))
                .findFirst().orElse(null);
    }

    private boolean checkGroupsAreOrdered()
    {
        List<PhasedVariantGroup> groups = getPhasedGroups();

        for(int i = 0; i < groups.size() - 1; ++i)
        {
            PhasedVariantGroup current = groups.get(i);
            PhasedVariantGroup next = groups.get(i + 1);

            if(next.posVariantMin() < current.posVariantMin())
                return false;
        }

        return true;
    }

    private SageVariant createVariant(int position)
    {
        VariantHotspot variant = ImmutableVariantHotspotImpl.builder()
                .chromosome("1")
                .position(position)
                .ref("A")
                .alt("C").build();

        List<ReadContextCounter> tumorCounters = Lists.newArrayList(createReadCounter(position));

        Candidate candidate = new Candidate(
                VariantTier.HIGH_CONFIDENCE, variant, tumorCounters.get(0).readContext(),
                100, 1, 1);

        List<ReadContextCounter> normalCounters = Lists.newArrayList();
        return new SageVariant(candidate, normalCounters, tumorCounters);
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
