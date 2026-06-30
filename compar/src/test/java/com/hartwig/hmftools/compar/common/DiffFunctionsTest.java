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
    private static final CategoryType CATEGORY = CategoryType.CDR3_LOCUS_SUMMARY;

    @Test
    public void checkDiffIntWithAbsoluteThreshold()
    {
        FieldConfig fieldConfig = new FieldConfig();

        fieldConfig.addFieldThreshold(CATEGORY, FIELD, 10, -1);

        assertNoDiff(0, 0, fieldConfig);
        assertNoDiff(100, 97, fieldConfig);
        assertNoDiff(-10, -4, fieldConfig);

        assertDiff(-10, 10, fieldConfig);
        assertDiff(10, -10, fieldConfig);
        assertDiff(1, 100, fieldConfig);
        assertDiff(100, 1, fieldConfig);

        assertNoDiff(0, 10, fieldConfig);
        assertNoDiff(10, 0, fieldConfig);
        assertDiff(0, 11, fieldConfig);
        assertDiff(11, 0, fieldConfig);

        assertNoDiff(90, 100, fieldConfig);
        assertNoDiff(100, 90, fieldConfig);
        assertDiff(89, 100, fieldConfig);
        assertDiff(100, 89, fieldConfig);

        assertNoDiff(-200, -190, fieldConfig);
        assertNoDiff(-190, -200, fieldConfig);
        assertDiff(-200, -189, fieldConfig);
        assertDiff(-189, -200, fieldConfig);

        assertNoDiff(-5, 4, fieldConfig);
        assertNoDiff(4, -5, fieldConfig);
        assertDiff(-5, 6, fieldConfig);
        assertDiff(6, -5, fieldConfig);
    }

    @Test
    public void checkDiffIntWithRelativeThreshold()
    {
        FieldConfig fieldConfig = new FieldConfig();

        fieldConfig.addFieldThreshold(CATEGORY, FIELD, -1, 0.2);

        assertNoDiff(0, 0, fieldConfig);
        assertNoDiff(100, 97, fieldConfig);
        assertNoDiff(-10, -9, fieldConfig);

        assertDiff(-1, 1, fieldConfig);
        assertDiff(1, -1, fieldConfig);

        assertDiff(0, 1, fieldConfig);
        assertDiff(1, 0, fieldConfig);
        assertDiff(0, -1, fieldConfig);
        assertDiff(-1, 0, fieldConfig);

        assertNoDiff(81, 100, fieldConfig);
        assertNoDiff(100, 81, fieldConfig);
        assertDiff(79, 100, fieldConfig);
        assertDiff(100, 79, fieldConfig);

        assertNoDiff(-1000, -810, fieldConfig);
        assertNoDiff(-810, -1000, fieldConfig);
        assertDiff(-1000, -790, fieldConfig);
        assertDiff(-790, -1000, fieldConfig);
    }

    @Test
    public void checkDiffIntWithEmptyThreshold()
    {
        FieldConfig fieldConfig = new FieldConfig();

        assertNoDiff(0, 0, fieldConfig);
        assertNoDiff(100, 100, fieldConfig);
        assertNoDiff(-10, -10, fieldConfig);

        assertDiff(0, 1, fieldConfig);
        assertDiff(1, 0, fieldConfig);

        assertDiff(99, 100, fieldConfig);
        assertDiff(100, 99, fieldConfig);

        assertDiff(-100, -99, fieldConfig);
        assertDiff(-99, -100, fieldConfig);

        assertDiff(-10, 10, fieldConfig);
        assertDiff(10, -10, fieldConfig);
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
        FieldConfig fieldConfig = new FieldConfig();

        assertNoDiff(0., 0., fieldConfig);
        assertNoDiff(100., 100., fieldConfig);
        assertNoDiff(-10., -10., fieldConfig);

        assertDiff(0., 1.1, fieldConfig);
        assertDiff(1.1, 0., fieldConfig);
        assertNoDiff(0., 0.9, fieldConfig);
        assertNoDiff(0.9, 0., fieldConfig);

        assertDiff(-0., -1.1, fieldConfig);
        assertDiff(-1.1, -0., fieldConfig);
        assertNoDiff(-0., -0.9, fieldConfig);
        assertNoDiff(-0.9, -0., fieldConfig);

        assertDiff(89., 100., fieldConfig);
        assertDiff(100., 89., fieldConfig);
        assertNoDiff(91., 100., fieldConfig);
        assertNoDiff(100., 91., fieldConfig);

        assertDiff(-100., -89., fieldConfig);
        assertDiff(-89., -100., fieldConfig);
        assertNoDiff(-100., -91., fieldConfig);
        assertNoDiff(-91., -100., fieldConfig);

        assertDiff(-10., 10., fieldConfig);
        assertDiff(10., -10., fieldConfig);
    }

    @Test
    public void checkDiffDoubleWithThreshold()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.addFieldThreshold(CATEGORY, FIELD, 20, -1);

        assertNoDiff(0., 0., fieldConfig);
        assertNoDiff(100., 100., fieldConfig);
        assertNoDiff(-10., -10., fieldConfig);

        assertDiff(0., 20.5, fieldConfig);
        assertDiff(20.5, 0., fieldConfig);
        assertNoDiff(0., 19.5, fieldConfig);
        assertNoDiff(19.5, 0., fieldConfig);

        assertDiff(-0., -20.5, fieldConfig);
        assertDiff(-20.5, -0., fieldConfig);
        assertNoDiff(-0., -19.5, fieldConfig);
        assertNoDiff(-19.5, -0., fieldConfig);

        assertDiff(79., 100., fieldConfig);
        assertDiff(100., 79., fieldConfig);
        assertNoDiff(81., 100., fieldConfig);
        assertNoDiff(100., 81., fieldConfig);

        assertDiff(-100., -79., fieldConfig);
        assertDiff(-79., -100., fieldConfig);
        assertNoDiff(-100., -81., fieldConfig);
        assertNoDiff(-81., -100., fieldConfig);

        assertDiff(-10., 10.5, fieldConfig);
        assertDiff(10.5, -10., fieldConfig);
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

    private static void assertNoDiff(final double refValue, final double newValue, final FieldConfig fieldConfig)
    {
        List<String> diffs = new ArrayList<>();
        assertFalse(checkDiff(diffs, FIELD, refValue, newValue, CATEGORY, fieldConfig));
        assertTrue(diffs.isEmpty());
    }

    private static void assertNoDiff(final int refValue, final int newValue, final FieldConfig fieldConfig)
    {
        List<String> diffs = new ArrayList<>();
        assertFalse(checkDiff(diffs, FIELD, refValue, newValue, CATEGORY, fieldConfig));
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

    private static void assertDiff(final double refValue, final double newValue, final FieldConfig fieldConfig)
    {
        List<String> diffs = new ArrayList<>();
        assertTrue(checkDiff(diffs, FIELD, refValue, newValue, CATEGORY, fieldConfig));
        assertEquals(1, diffs.size());
        assertEquals(FIELD, extractFieldNameFromDifference(diffs.get(0)));
    }

    private static void assertDiff(final int refValue, final int newValue, final FieldConfig fieldConfig)
    {
        List<String> diffs = new ArrayList<>();
        assertTrue(checkDiff(diffs, FIELD, refValue, newValue, CATEGORY, fieldConfig));
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
