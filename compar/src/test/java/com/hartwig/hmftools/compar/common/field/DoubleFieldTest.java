package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.InvalidDataItem;

import org.junit.Test;

public class DoubleFieldTest
{
    private static final String FIELD_NAME = "DoubleField";

    private static DoubleField field(final Double absoluteThreshold, final Double percentThreshold)
    {
        return new DoubleField(FIELD_NAME, i -> ((TestFieldItem<Double>) i).Value, true, absoluteThreshold,
                percentThreshold, "%.1f");
    }

    private static boolean hasDiff(final DoubleField field, final double oldValue, final double newValue)
    {
        return field.hasDiff(new TestFieldItem<>(oldValue), new TestFieldItem<>(newValue));
    }

    @Test
    public void nameIsComparedAndThresholdsReflectConstructorArgs()
    {
        DoubleField field = field(1., 0.2);
        assertEquals(FIELD_NAME, field.name());
        assertTrue(field.isCompared());
        assertEquals(1., field.absoluteThreshold(), 0.0);
        assertEquals(0.2, field.percentThreshold(), 0.0);
    }

    @Test
    public void thresholdsAreNullByDefault()
    {
        DoubleField field = field(null, null);
        assertNull(field.absoluteThreshold());
        assertNull(field.percentThreshold());
    }

    @Test
    public void displayValueUsesFormatString()
    {
        DoubleField field = field(null, null);
        assertEquals("5.0", field.displayValue(new TestFieldItem<>(5.0)));
        assertEquals("-3.1", field.displayValue(new TestFieldItem<>(-3.14)));
    }

    @Test
    public void displayValueIsEmptyForInvalidItem()
    {
        DoubleField field = field(null, null);
        assertEquals("", field.displayValue(new InvalidDataItem(CategoryType.PURITY)));
    }

    @Test
    public void equalValuesNeverDiffRegardlessOfThresholds()
    {
        DoubleField field = field(null, null);
        assertFalse(hasDiff(field, 5.0, 5.0));
        assertFalse(hasDiff(field, 0.0, 0.0));
        assertFalse(hasDiff(field, -5.0, -5.0));
        assertFalse(hasDiff(field, -0.0, 0.0));
    }

    @Test
    public void withoutThresholdsAnyDifferenceIsADiff()
    {
        DoubleField field = field(null, null);
        assertTrue(hasDiff(field, 5.0, 5.01));
        assertTrue(hasDiff(field, -5.0, -5.01));
        assertTrue(hasDiff(field, -5.0, 5.0));
    }

    @Test
    public void absoluteThresholdOnlyRequiresAbsoluteDiffToExceedThreshold()
    {
        DoubleField field = field(1., null);
        assertFalse(hasDiff(field, 10.0, 10.9));
        assertTrue(hasDiff(field, 10.0, 11.1));
    }

    @Test
    public void absoluteThresholdOnlyWorksWithNegativeAndMixedSignValues()
    {
        DoubleField field = field(1., null);
        assertFalse(hasDiff(field, -10.0, -10.9));
        assertTrue(hasDiff(field, -10.0, -11.1));

        // absDiff = 10 > 1 -> diff, even though the values have opposite signs
        assertTrue(hasDiff(field, -5.0, 5.0));
    }

    @Test
    public void percentThresholdOnlyRequiresRelativeDiffToExceedThreshold()
    {
        DoubleField field = field(null, 0.2);
        assertFalse(hasDiff(field, 81.0, 100.0));
        assertTrue(hasDiff(field, 79.0, 100.0));
    }

    @Test
    public void percentThresholdOnlyWorksWithNegativeAndMixedSignValues()
    {
        DoubleField field = field(null, 0.2);
        assertFalse(hasDiff(field, -81.0, -100.0));
        assertTrue(hasDiff(field, -79.0, -100.0));

        // absDiff = 100, relDiff = 100/50 = 2.0 -> diff, even though the values have opposite signs
        assertTrue(hasDiff(field, -50.0, 50.0));
    }

    @Test
    public void bothThresholdsRequireBothToBeExceeded()
    {
        DoubleField field = field(1., 0.2);

        // absolute and percent both exceeded -> diff
        assertTrue(hasDiff(field, 10.0, 100.0));

        // absolute exceeded, percent not -> no diff
        assertFalse(hasDiff(field, 1000.0, 1011.0));

        // percent exceeded, absolute not -> no diff
        assertFalse(hasDiff(field, 1.0, 1.3));

        // neither exceeded -> no diff
        assertFalse(hasDiff(field, 5.0, 5.5));
    }

    @Test
    public void bothThresholdsRequireBothToBeExceededWithNegativeAndMixedSignValues()
    {
        DoubleField field = field(1., 0.2);

        // absolute and percent both exceeded, negative values -> diff
        assertTrue(hasDiff(field, -10.0, -100.0));

        // absolute exceeded, percent not, negative values -> no diff
        assertFalse(hasDiff(field, -1000.0, -1011.0));

        // percent exceeded, absolute not, negative values -> no diff
        assertFalse(hasDiff(field, -1.0, -1.3));

        // absolute and percent both exceeded, opposite signs -> diff
        assertTrue(hasDiff(field, -50.0, 50.0));
    }

    @Test
    public void withComparedUpdatesIsComparedAndPreservesThresholds()
    {
        DoubleField field = field(1., 0.2);
        DoubleField updated = (DoubleField) field.withCompared(false);

        assertFalse(updated.isCompared());
        assertEquals(1., updated.absoluteThreshold(), 0.0);
        assertEquals(0.2, updated.percentThreshold(), 0.0);
    }

    @Test
    public void withAbsoluteThresholdUpdatesOnlyAbsoluteThreshold()
    {
        DoubleField field = field(1., 0.2);
        DoubleField updated = (DoubleField) field.withAbsoluteThreshold(5.);

        assertTrue(updated.isCompared());
        assertEquals(5., updated.absoluteThreshold(), 0.0);
        assertEquals(0.2, updated.percentThreshold(), 0.0);
    }

    @Test
    public void withPercentThresholdUpdatesOnlyPercentThreshold()
    {
        DoubleField field = field(1., 0.2);
        DoubleField updated = (DoubleField) field.withPercentThreshold(0.5);

        assertTrue(updated.isCompared());
        assertEquals(1., updated.absoluteThreshold(), 0.0);
        assertEquals(0.5, updated.percentThreshold(), 0.0);
    }
}
