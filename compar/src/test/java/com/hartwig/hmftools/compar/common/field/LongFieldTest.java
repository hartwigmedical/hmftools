package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class LongFieldTest
{
    private static final String FIELD_NAME = "LongField";

    private static LongField field(final Double absoluteThreshold, final Double percentThreshold)
    {
        return new LongField(FIELD_NAME, i -> ((TestFieldItem<Long>) i).Value, true,
                absoluteThreshold, percentThreshold);
    }

    private static boolean hasDiff(final LongField field, final long oldValue, final long newValue)
    {
        return field.hasDiff(new TestFieldItem<>(oldValue), new TestFieldItem<>(newValue));
    }

    @Test
    public void nameIsComparedAndThresholdsReflectConstructorArgs()
    {
        LongField field = field(10., 0.2);
        assertEquals(FIELD_NAME, field.name());
        assertTrue(field.isCompared());
        assertEquals(10., field.absoluteThreshold(), 0.0);
        assertEquals(0.2, field.percentThreshold(), 0.0);
    }

    @Test
    public void thresholdsAreNullByDefault()
    {
        LongField field = field(null, null);
        assertNull(field.absoluteThreshold());
        assertNull(field.percentThreshold());
    }

    @Test
    public void displayValueIsLongToString()
    {
        LongField field = field(null, null);
        assertEquals("5000000000", field.displayValue(new TestFieldItem<>(5_000_000_000L)));
        assertEquals("-3", field.displayValue(new TestFieldItem<>(-3L)));
    }

    @Test
    public void equalValuesNeverDiffRegardlessOfThresholds()
    {
        LongField field = field(null, null);
        assertFalse(hasDiff(field, 5L, 5L));
        assertFalse(hasDiff(field, 0L, 0L));
        assertFalse(hasDiff(field, -5L, -5L));
    }

    @Test
    public void withoutThresholdsAnyDifferenceIsADiff()
    {
        LongField field = field(null, null);
        assertTrue(hasDiff(field, 5L, 6L));
        assertTrue(hasDiff(field, 100L, 99L));
        assertTrue(hasDiff(field, -5L, -6L));
        assertTrue(hasDiff(field, -5L, 5L));
    }

    @Test
    public void absoluteThresholdOnlyRequiresAbsoluteDiffToExceedThreshold()
    {
        LongField field = field(10., null);
        assertFalse(hasDiff(field, 0L, 10L));
        assertTrue(hasDiff(field, 0L, 11L));
        assertTrue(hasDiff(field, -10L, 10L));
    }

    @Test
    public void absoluteThresholdOnlyWorksWithNegativeValues()
    {
        LongField field = field(10., null);
        assertFalse(hasDiff(field, 0L, -10L));
        assertTrue(hasDiff(field, 0L, -11L));
        assertFalse(hasDiff(field, -100L, -90L));
        assertTrue(hasDiff(field, -100L, -89L));
    }

    @Test
    public void percentThresholdOnlyRequiresRelativeDiffToExceedThreshold()
    {
        LongField field = field(null, 0.2);
        assertFalse(hasDiff(field, 81L, 100L));
        assertTrue(hasDiff(field, 79L, 100L));
    }

    @Test
    public void percentThresholdOnlyWorksWithNegativeAndMixedSignValues()
    {
        LongField field = field(null, 0.2);
        assertFalse(hasDiff(field, -81L, -100L));
        assertTrue(hasDiff(field, -79L, -100L));

        // absDiff = 100, relDiff = 100/50 = 2.0 -> diff, even though the values have opposite signs
        assertTrue(hasDiff(field, -50L, 50L));
    }

    @Test
    public void bothThresholdsRequireBothToBeExceeded()
    {
        LongField field = field(10., 0.2);

        // absolute and percent both exceeded -> diff
        assertTrue(hasDiff(field, 10L, 100L));

        // absolute exceeded, percent not -> no diff
        assertFalse(hasDiff(field, 1000L, 1011L));

        // percent exceeded, absolute not -> no diff
        assertFalse(hasDiff(field, 1L, 2L));

        // neither exceeded -> no diff
        assertFalse(hasDiff(field, 5L, 6L));
    }

    @Test
    public void bothThresholdsRequireBothToBeExceededWithNegativeAndMixedSignValues()
    {
        LongField field = field(10., 0.2);

        // absolute and percent both exceeded, negative values -> diff
        assertTrue(hasDiff(field, -10L, -100L));

        // absolute exceeded, percent not, negative values -> no diff
        assertFalse(hasDiff(field, -1000L, -1011L));

        // percent exceeded, absolute not, negative values -> no diff
        assertFalse(hasDiff(field, -1L, -2L));

        // absolute and percent both exceeded, opposite signs -> diff
        assertTrue(hasDiff(field, -50L, 50L));
    }

    @Test
    public void withComparedUpdatesIsComparedAndPreservesThresholds()
    {
        LongField field = field(10., 0.2);
        LongField updated = (LongField) field.withCompared(false);

        assertFalse(updated.isCompared());
        assertEquals(10., updated.absoluteThreshold(), 0.0);
        assertEquals(0.2, updated.percentThreshold(), 0.0);
    }

    @Test
    public void withAbsoluteThresholdUpdatesOnlyAbsoluteThreshold()
    {
        LongField field = field(10., 0.2);
        LongField updated = (LongField) field.withAbsoluteThreshold(5.);

        assertTrue(updated.isCompared());
        assertEquals(5., updated.absoluteThreshold(), 0.0);
        assertEquals(0.2, updated.percentThreshold(), 0.0);
    }

    @Test
    public void withPercentThresholdUpdatesOnlyPercentThreshold()
    {
        LongField field = field(10., 0.2);
        LongField updated = (LongField) field.withPercentThreshold(0.5);

        assertTrue(updated.isCompared());
        assertEquals(10., updated.absoluteThreshold(), 0.0);
        assertEquals(0.5, updated.percentThreshold(), 0.0);
    }
}
