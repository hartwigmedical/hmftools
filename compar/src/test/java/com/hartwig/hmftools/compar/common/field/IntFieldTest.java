package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class IntFieldTest
{
    private static final String FIELD_NAME = "IntField";

    private static IntField field(final Double absoluteThreshold, final Double percentThreshold)
    {
        return new IntField(FIELD_NAME, i -> ((TestFieldItem<Integer>) i).Value, true, absoluteThreshold, percentThreshold);
    }

    private static boolean hasDiff(final IntField field, final int oldValue, final int newValue)
    {
        return field.hasDiff(new TestFieldItem<>(oldValue), new TestFieldItem<>(newValue));
    }

    @Test
    public void nameIsComparedAndThresholdsReflectConstructorArgs()
    {
        IntField field = field(10., 0.2);
        assertEquals(FIELD_NAME, field.name());
        assertTrue(field.isCompared());
        assertEquals(10., field.absoluteThreshold(), 0.0);
        assertEquals(0.2, field.percentThreshold(), 0.0);
    }

    @Test
    public void thresholdsAreNullByDefault()
    {
        IntField field = field(null, null);
        assertNull(field.absoluteThreshold());
        assertNull(field.percentThreshold());
    }

    @Test
    public void displayValueIsIntegerToString()
    {
        IntField field = field(null, null);
        assertEquals("5", field.displayValue(new TestFieldItem<>(5)));
        assertEquals("-3", field.displayValue(new TestFieldItem<>(-3)));
    }

    @Test
    public void equalValuesNeverDiffRegardlessOfThresholds()
    {
        IntField field = field(null, null);
        assertFalse(hasDiff(field, 5, 5));
        assertFalse(hasDiff(field, 0, 0));
        assertFalse(hasDiff(field, -5, -5));
    }

    @Test
    public void withoutThresholdsAnyDifferenceIsADiff()
    {
        IntField field = field(null, null);
        assertTrue(hasDiff(field, 5, 6));
        assertTrue(hasDiff(field, 100, 99));
        assertTrue(hasDiff(field, -5, -6));
        assertTrue(hasDiff(field, -5, 5));
    }

    @Test
    public void absoluteThresholdOnlyRequiresAbsoluteDiffToExceedThreshold()
    {
        IntField field = field(10., null);
        assertFalse(hasDiff(field, 0, 10));
        assertTrue(hasDiff(field, 0, 11));
        assertTrue(hasDiff(field, -10, 10));
    }

    @Test
    public void absoluteThresholdOnlyWorksWithNegativeValues()
    {
        IntField field = field(10., null);
        assertFalse(hasDiff(field, 0, -10));
        assertTrue(hasDiff(field, 0, -11));
        assertFalse(hasDiff(field, -100, -90));
        assertTrue(hasDiff(field, -100, -89));
    }

    @Test
    public void percentThresholdOnlyRequiresRelativeDiffToExceedThreshold()
    {
        IntField field = field(null, 0.2);
        assertFalse(hasDiff(field, 81, 100));
        assertTrue(hasDiff(field, 79, 100));
    }

    @Test
    public void percentThresholdOnlyWorksWithNegativeAndMixedSignValues()
    {
        IntField field = field(null, 0.2);
        assertFalse(hasDiff(field, -81, -100));
        assertTrue(hasDiff(field, -79, -100));

        // absDiff = 100, relDiff = 100/50 = 2.0 -> diff, even though the values have opposite signs
        assertTrue(hasDiff(field, -50, 50));
    }

    @Test
    public void bothThresholdsRequireBothToBeExceeded()
    {
        IntField field = field(10., 0.2);

        // absolute and percent both exceeded -> diff
        assertTrue(hasDiff(field, 10, 100));

        // absolute exceeded, percent not -> no diff
        assertFalse(hasDiff(field, 1000, 1011));

        // percent exceeded, absolute not -> no diff
        assertFalse(hasDiff(field, 1, 2));

        // neither exceeded -> no diff
        assertFalse(hasDiff(field, 5, 6));
    }

    @Test
    public void bothThresholdsRequireBothToBeExceededWithNegativeAndMixedSignValues()
    {
        IntField field = field(10., 0.2);

        // absolute and percent both exceeded, negative values -> diff
        assertTrue(hasDiff(field, -10, -100));

        // absolute exceeded, percent not, negative values -> no diff
        assertFalse(hasDiff(field, -1000, -1011));

        // percent exceeded, absolute not, negative values -> no diff
        assertFalse(hasDiff(field, -1, -2));

        // absolute and percent both exceeded, opposite signs -> diff
        assertTrue(hasDiff(field, -50, 50));
    }

    @Test
    public void withComparedUpdatesIsComparedAndPreservesThresholds()
    {
        IntField field = field(10., 0.2);
        IntField updated = (IntField) field.withCompared(false);

        assertFalse(updated.isCompared());
        assertEquals(10., updated.absoluteThreshold(), 0.0);
        assertEquals(0.2, updated.percentThreshold(), 0.0);
    }

    @Test
    public void withAbsoluteThresholdUpdatesOnlyAbsoluteThreshold()
    {
        IntField field = field(10., 0.2);
        IntField updated = (IntField) field.withAbsoluteThreshold(5.);

        assertTrue(updated.isCompared());
        assertEquals(5., updated.absoluteThreshold(), 0.0);
        assertEquals(0.2, updated.percentThreshold(), 0.0);
    }

    @Test
    public void withPercentThresholdUpdatesOnlyPercentThreshold()
    {
        IntField field = field(10., 0.2);
        IntField updated = (IntField) field.withPercentThreshold(0.5);

        assertTrue(updated.isCompared());
        assertEquals(10., updated.absoluteThreshold(), 0.0);
        assertEquals(0.5, updated.percentThreshold(), 0.0);
    }
}
