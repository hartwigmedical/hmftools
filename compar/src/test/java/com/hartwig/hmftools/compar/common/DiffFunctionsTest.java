package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.ComparTestUtil.extractFieldNameFromDifference;
import static com.hartwig.hmftools.compar.common.DiffFunctions.FILTER_DIFF;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkFilterDiffs;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

public class DiffFunctionsTest
{

    public static final String FIELD = "Field";

    @Test
    public void checkDiffIntWithAbsoluteThreshold()
    {
        DiffThresholds diffThresholds = new DiffThresholds();

        diffThresholds.addFieldThreshold(FIELD, 10, -1);

        assertNoDiff(0, 0, diffThresholds);
        assertNoDiff(100, 97, diffThresholds);
        assertNoDiff(-10, -4, diffThresholds);

        assertDiff(-10, 10, diffThresholds);
        assertDiff(10, -10, diffThresholds);
        assertDiff(1, 100, diffThresholds);
        assertDiff(100, 1, diffThresholds);

        assertNoDiff(0, 10, diffThresholds);
        assertNoDiff(10, 0, diffThresholds);
        assertDiff(0, 11, diffThresholds);
        assertDiff(11, 0, diffThresholds);

        assertNoDiff(90, 100, diffThresholds);
        assertNoDiff(100, 90, diffThresholds);
        assertDiff(89, 100, diffThresholds);
        assertDiff(100, 89, diffThresholds);

        assertNoDiff(-200, -190, diffThresholds);
        assertNoDiff(-190, -200, diffThresholds);
        assertDiff(-200, -189, diffThresholds);
        assertDiff(-189, -200, diffThresholds);

        assertNoDiff(-5, 4, diffThresholds);
        assertNoDiff(4, -5, diffThresholds);
        assertDiff(-5, 6, diffThresholds);
        assertDiff(6, -5, diffThresholds);
    }

    @Test
    public void checkDiffIntWithRelativeThreshold()
    {
        DiffThresholds diffThresholds = new DiffThresholds();

        diffThresholds.addFieldThreshold(FIELD, -1, 0.2);

        assertNoDiff(0, 0, diffThresholds);
        assertNoDiff(100, 97, diffThresholds);
        assertNoDiff(-10, -9, diffThresholds);

        assertDiff(-1, 1, diffThresholds);
        assertDiff(1, -1, diffThresholds);

        assertDiff(0, 1, diffThresholds);
        assertDiff(1, 0, diffThresholds);
        assertDiff(0, -1, diffThresholds);
        assertDiff(-1, 0, diffThresholds);

        assertNoDiff(81, 100, diffThresholds);
        assertNoDiff(100, 81, diffThresholds);
        assertDiff(79, 100, diffThresholds);
        assertDiff(100, 79, diffThresholds);

        assertNoDiff(-1000, -810, diffThresholds);
        assertNoDiff(-810, -1000, diffThresholds);
        assertDiff(-1000, -790, diffThresholds);
        assertDiff(-790, -1000, diffThresholds);
    }

    @Test
    public void checkDiffIntWithEmptyThreshold()
    {
        DiffThresholds diffThresholds = new DiffThresholds();

        assertNoDiff(0, 0, diffThresholds);
        assertNoDiff(100, 100, diffThresholds);
        assertNoDiff(-10, -10, diffThresholds);

        assertDiff(0, 1, diffThresholds);
        assertDiff(1, 0, diffThresholds);

        assertDiff(99, 100, diffThresholds);
        assertDiff(100, 99, diffThresholds);

        assertDiff(-100, -99, diffThresholds);
        assertDiff(-99, -100, diffThresholds);

        assertDiff(-10, 10, diffThresholds);
        assertDiff(10, -10, diffThresholds);
    }

    @Test
    public void checkDiffIntWithoutThreshold()
    {
        assertNoDiff(0, 0);
        assertNoDiff(100, 100);
        assertNoDiff(-10, -10);

        assertDiff(0, 1);
        assertDiff(1, 0);

        assertDiff(99, 100);
        assertDiff(100, 99);

        assertDiff(-100, -99);
        assertDiff(-99, -100);

        assertDiff(-10, 10);
        assertDiff(10, -10);
    }

    @Test
    public void checkDiffDoubleWithEmptyThreshold()
    {
        DiffThresholds diffThresholds = new DiffThresholds();

        assertNoDiff(0., 0., diffThresholds);
        assertNoDiff(100., 100., diffThresholds);
        assertNoDiff(-10., -10., diffThresholds);

        assertDiff(0., 1.1, diffThresholds);
        assertDiff(1.1, 0., diffThresholds);
        assertNoDiff(0., 0.9, diffThresholds);
        assertNoDiff(0.9, 0., diffThresholds);

        assertDiff(-0., -1.1, diffThresholds);
        assertDiff(-1.1, -0., diffThresholds);
        assertNoDiff(-0., -0.9, diffThresholds);
        assertNoDiff(-0.9, -0., diffThresholds);

        assertDiff(89., 100., diffThresholds);
        assertDiff(100., 89., diffThresholds);
        assertNoDiff(91., 100., diffThresholds);
        assertNoDiff(100., 91., diffThresholds);

        assertDiff(-100., -89., diffThresholds);
        assertDiff(-89., -100., diffThresholds);
        assertNoDiff(-100., -91., diffThresholds);
        assertNoDiff(-91., -100., diffThresholds);

        assertDiff(-10., 10., diffThresholds);
        assertDiff(10., -10., diffThresholds);
    }

    @Test
    public void checkDiffDoubleWithThreshold()
    {
        DiffThresholds diffThresholds = new DiffThresholds();
        diffThresholds.addFieldThreshold(FIELD, 20, -1);

        assertNoDiff(0., 0., diffThresholds);
        assertNoDiff(100., 100., diffThresholds);
        assertNoDiff(-10., -10., diffThresholds);

        assertDiff(0., 20.5, diffThresholds);
        assertDiff(20.5, 0., diffThresholds);
        assertNoDiff(0., 19.5, diffThresholds);
        assertNoDiff(19.5, 0., diffThresholds);

        assertDiff(-0., -20.5, diffThresholds);
        assertDiff(-20.5, -0., diffThresholds);
        assertNoDiff(-0., -19.5, diffThresholds);
        assertNoDiff(-19.5, -0., diffThresholds);

        assertDiff(79., 100., diffThresholds);
        assertDiff(100., 79., diffThresholds);
        assertNoDiff(81., 100., diffThresholds);
        assertNoDiff(100., 81., diffThresholds);

        assertDiff(-100., -79., diffThresholds);
        assertDiff(-79., -100., diffThresholds);
        assertNoDiff(-100., -81., diffThresholds);
        assertNoDiff(-81., -100., diffThresholds);

        assertDiff(-10., 10.5, diffThresholds);
        assertDiff(10.5, -10., diffThresholds);
    }

    @Test
    public void checkDiffBoolean()
    {
        assertNoDiff(true, true);
        assertNoDiff(false, false);

        assertDiff(false, true);
        assertDiff(true, false);
    }

    @Test
    public void checkDiffString()
    {
        assertNoDiff("", "");
        assertNoDiff("TEST", "TEST");

        assertDiff("", "TEST");
        assertDiff("TEST", "");

        assertDiff("test", "TEST");
        assertDiff("TEST", "test");

        assertDiff("ABC", "TEST");
    }

    @Test
    public void checkDiffFilters()
    {
        assertNoDiff(new HashSet<>(), new HashSet<>());
        assertNoDiff(Set.of("HI"), Set.of("HI"));

        assertDiff(Set.of("HI"), new HashSet<>());
        assertDiff(new HashSet<>(), Set.of("HI"));

        assertDiff(Set.of("HI"), Set.of("hi"));
        assertDiff(Set.of("hi"), Set.of("HI"));

        assertDiff(Set.of("hello", "world"), Set.of("hello"));
        assertDiff(Set.of("hello"), Set.of("hello", "world"));

        assertDiff(Set.of("hi", "world"), Set.of("hello", "world"));
    }

    private static void assertNoDiff(final double refValue, final double newValue, final DiffThresholds diffThresholds)
    {
        List<String> diffs = new ArrayList<>();
        assertFalse(checkDiff(diffs, FIELD, refValue, newValue, diffThresholds));
        assertTrue(diffs.isEmpty());
    }

    private static void assertNoDiff(final int refValue, final int newValue, final DiffThresholds diffThresholds)
    {
        List<String> diffs = new ArrayList<>();
        assertFalse(checkDiff(diffs, FIELD, refValue, newValue, diffThresholds));
        assertTrue(diffs.isEmpty());
    }

    private static void assertNoDiff(final int refValue, final int newValue)
    {
        List<String> diffs = new ArrayList<>();
        assertFalse(checkDiff(diffs, FIELD, refValue, newValue));
        assertTrue(diffs.isEmpty());
    }

    private static void assertNoDiff(final boolean refValue, final boolean newValue)
    {
        List<String> diffs = new ArrayList<>();
        assertFalse(checkDiff(diffs, FIELD, refValue, newValue));
        assertTrue(diffs.isEmpty());
    }

    private static void assertNoDiff(final String refValue, final String newValue)
    {
        List<String> diffs = new ArrayList<>();
        assertFalse(checkDiff(diffs, FIELD, refValue, newValue));
        assertTrue(diffs.isEmpty());
    }

    private static void assertNoDiff(final Set<String> refValue, final Set<String> newValue)
    {
        List<String> diffs = new ArrayList<>();
        checkFilterDiffs(refValue, newValue, diffs);
        assertTrue(diffs.isEmpty());
    }

    private static void assertDiff(final double refValue, final double newValue, final DiffThresholds diffThresholds)
    {
        List<String> diffs = new ArrayList<>();
        assertTrue(checkDiff(diffs, FIELD, refValue, newValue, diffThresholds));
        assertEquals(1, diffs.size());
        assertEquals(FIELD, extractFieldNameFromDifference(diffs.get(0)));
    }

    private static void assertDiff(final int refValue, final int newValue, final DiffThresholds diffThresholds)
    {
        List<String> diffs = new ArrayList<>();
        assertTrue(checkDiff(diffs, FIELD, refValue, newValue, diffThresholds));
        assertEquals(1, diffs.size());
        assertEquals(FIELD, extractFieldNameFromDifference(diffs.get(0)));
    }

    private static void assertDiff(final int refValue, final int newValue)
    {
        List<String> diffs = new ArrayList<>();
        assertTrue(checkDiff(diffs, FIELD, refValue, newValue));
        assertEquals(1, diffs.size());
        assertEquals(FIELD, extractFieldNameFromDifference(diffs.get(0)));
    }

    private static void assertDiff(final boolean refValue, final boolean newValue)
    {
        List<String> diffs = new ArrayList<>();
        assertTrue(checkDiff(diffs, FIELD, refValue, newValue));
        assertEquals(1, diffs.size());
        assertEquals(FIELD, extractFieldNameFromDifference(diffs.get(0)));
    }

    private static void assertDiff(final String refValue, final String newValue)
    {
        List<String> diffs = new ArrayList<>();
        assertTrue(checkDiff(diffs, FIELD, refValue, newValue));
        assertEquals(1, diffs.size());
        assertEquals(FIELD, extractFieldNameFromDifference(diffs.get(0)));
    }

    private static void assertDiff(final Set<String> refValue, final Set<String> newValue)
    {
        List<String> diffs = new ArrayList<>();
        checkFilterDiffs(refValue, newValue, diffs);
        assertEquals(1, diffs.size());
        assertEquals(FILTER_DIFF, extractFieldNameFromDifference(diffs.get(0)));
    }
}
