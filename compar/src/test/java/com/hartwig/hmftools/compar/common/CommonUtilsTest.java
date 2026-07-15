package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonChromosome;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.common.MatchLevel.DETAILED;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.common.MismatchType.FULL_MATCH;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.OLD_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;
import static com.hartwig.hmftools.compar.common.SourceType.NEW;
import static com.hartwig.hmftools.compar.common.SourceType.OLD;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.ComparableItem;

import org.junit.Test;

public class CommonUtilsTest
{
    private static class TestComparableItem implements ComparableItem
    {
        private final boolean mReportable;
        private final boolean mIsPass;

        public TestComparableItem(final boolean reportable, final boolean isPass)
        {
            mReportable = reportable;
            mIsPass = isPass;
        }

        @Override
        public CategoryType category() { return CategoryType.DRIVER; }

        @Override
        public boolean matches(final ComparableItem other) { return true; }

        @Override
        public String key() { return "key"; }

        @Override
        public boolean reportable() { return mReportable; }

        @Override
        public boolean isPass() { return mIsPass; }
    }

    private static final List<String> NO_DIFFS = Lists.newArrayList();
    private static final List<String> SOME_DIFFS = Lists.newArrayList("field(old vs new)");

    @Test
    public void createMismatchFromDiffsIgnoresBothUncalledRegardlessOfDiffs()
    {
        TestComparableItem oldItem = new TestComparableItem(true, false);
        TestComparableItem newItem = new TestComparableItem(true, false);

        assertNull(createMismatchFromDiffs(oldItem, newItem, NO_DIFFS, DETAILED, false));
        assertNull(createMismatchFromDiffs(oldItem, newItem, NO_DIFFS, DETAILED, true));

        // even with diffs present, neither side counts as called so the difference is unimportant
        assertNull(createMismatchFromDiffs(oldItem, newItem, SOME_DIFFS, DETAILED, false));
        assertNull(createMismatchFromDiffs(oldItem, newItem, SOME_DIFFS, DETAILED, true));
    }

    @Test
    public void createMismatchFromDiffsIgnoresPerfectMatchUnlessIncludeMatches()
    {
        TestComparableItem oldItem = new TestComparableItem(true, true);
        TestComparableItem newItem = new TestComparableItem(true, true);

        assertNull(createMismatchFromDiffs(oldItem, newItem, NO_DIFFS, DETAILED, false));

        Mismatch mismatch = createMismatchFromDiffs(oldItem, newItem, NO_DIFFS, DETAILED, true);
        assertEquals(FULL_MATCH, mismatch.Type);
    }

    @Test
    public void createMismatchFromDiffsReturnsValueWhenBothCalledWithDiffs()
    {
        TestComparableItem oldItem = new TestComparableItem(true, true);
        TestComparableItem newItem = new TestComparableItem(true, true);

        Mismatch mismatch = createMismatchFromDiffs(oldItem, newItem, SOME_DIFFS, DETAILED, false);
        assertEquals(VALUE, mismatch.Type);

        // still VALUE when matches are also being included
        mismatch = createMismatchFromDiffs(oldItem, newItem, SOME_DIFFS, DETAILED, true);
        assertEquals(VALUE, mismatch.Type);
    }

    @Test
    public void createMismatchFromDiffsFlagsOldOnlyWhenOnlyOldCountsAsCalled()
    {
        TestComparableItem oldItem = new TestComparableItem(true, true);
        TestComparableItem newItem = new TestComparableItem(true, false);

        // this is the regression case: the two items were paired up despite one being an artificial
        // "not called" placeholder, and the always-compared fields happened to produce no diffs -
        // that must still surface as OLD_ONLY rather than being silently dropped
        Mismatch mismatch = createMismatchFromDiffs(oldItem, newItem, NO_DIFFS, DETAILED, false);
        assertEquals(OLD_ONLY, mismatch.Type);

        mismatch = createMismatchFromDiffs(oldItem, newItem, SOME_DIFFS, DETAILED, false);
        assertEquals(OLD_ONLY, mismatch.Type);
    }

    @Test
    public void createMismatchFromDiffsFlagsNewOnlyWhenOnlyNewCountsAsCalled()
    {
        TestComparableItem oldItem = new TestComparableItem(true, false);
        TestComparableItem newItem = new TestComparableItem(true, true);

        // regression case for the dropped NEW_ONLY AMP scenario: no diffs found, but the new side
        // is called and the old side is not, so this must still be reported as NEW_ONLY
        Mismatch mismatch = createMismatchFromDiffs(oldItem, newItem, NO_DIFFS, DETAILED, false);
        assertEquals(NEW_ONLY, mismatch.Type);

        mismatch = createMismatchFromDiffs(oldItem, newItem, SOME_DIFFS, DETAILED, false);
        assertEquals(NEW_ONLY, mismatch.Type);
    }

    @Test
    public void createMismatchFromDiffsUsesReportableAtReportableMatchLevel()
    {
        // at REPORTABLE match level, countsAsCalled() is based on reportable(), not isPass()
        TestComparableItem oldItem = new TestComparableItem(false, true);
        TestComparableItem newItem = new TestComparableItem(true, true);

        Mismatch mismatch = createMismatchFromDiffs(oldItem, newItem, NO_DIFFS, REPORTABLE, false);
        assertEquals(NEW_ONLY, mismatch.Type);

        // both reportable, no diffs -> ignored unless includeMatches
        TestComparableItem bothReportable1 = new TestComparableItem(true, false);
        TestComparableItem bothReportable2 = new TestComparableItem(true, false);
        assertNull(createMismatchFromDiffs(bothReportable1, bothReportable2, NO_DIFFS, REPORTABLE, false));
    }

    @Test
    public void emptyComparison()
    {
        FieldConfig fieldConfig = new FieldConfig();

        List<Mismatch> mismatches = new ArrayList<>();
        List<ComparableItem> refItems = new ArrayList<>();
        List<ComparableItem> newItems = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;
        boolean includeMatches = false;

        CommonUtils.compareItems(mismatches, matchLevel, fieldConfig, includeMatches, refItems, newItems);

        assertTrue(mismatches.isEmpty());
    }

    @Test
    public void testDetermineComparisonGenomePositionWithoutLiftover()
    {
        GenomeLiftoverCache lifoverCache = new GenomeLiftoverCache(true);

        assertEquals(new BasePosition("8", 10000),
                determineComparisonGenomePosition("8", 10000, OLD, false, lifoverCache));
        assertEquals(new BasePosition("chr10", 10000),
                determineComparisonGenomePosition("chr10", 10000, OLD, false, lifoverCache));
        assertEquals(new BasePosition("X", 10000),
                determineComparisonGenomePosition("X", 10000, NEW, false, lifoverCache));
        assertEquals(new BasePosition("chr3", 10000),
                determineComparisonGenomePosition("chr3", 10000, NEW, false, lifoverCache));
    }

    @Test
    public void testDetermineComparisonGenomePositionWithLiftover()
    {
        GenomeLiftoverCache lifoverCache = new GenomeLiftoverCache(true);

        assertTrue(lifoverCache.hasMappings());
        assertEquals(new BasePosition("chr8", 150000),
                determineComparisonGenomePosition("8", 100000, OLD, true, lifoverCache));
        assertEquals(new BasePosition("chr10", 100000),
                determineComparisonGenomePosition("chr10", 100000, OLD, true, lifoverCache));
        assertEquals(new BasePosition("X", 100000),
                determineComparisonGenomePosition("X", 100000, NEW, true, lifoverCache));
        assertEquals(new BasePosition("3", 141683),
                determineComparisonGenomePosition("chr3", 100000, NEW, true, lifoverCache));
    }

    @Test
    public void testDetermineComparisonChromosomeWithoutLiftover()
    {
        assertEquals("8", determineComparisonChromosome("8", false));
        assertEquals("chr10", determineComparisonChromosome("chr10", false));
        assertEquals("X", determineComparisonChromosome("X", false));
        assertEquals("chrY", determineComparisonChromosome("chrY", false));
  }

    @Test
    public void testDetermineComparisonChromosomeWithLiftover()
    {
        assertEquals("8", determineComparisonChromosome("8", true));
        assertEquals("10", determineComparisonChromosome("chr10", true));
        assertEquals("X", determineComparisonChromosome("X", true));
        assertEquals("Y", determineComparisonChromosome("chrY", true));
    }
}
