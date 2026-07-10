package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class BooleanFieldTest
{
    private static final String FIELD_NAME = "BoolField";

    private static BooleanField field(final boolean isCompared)
    {
        return new BooleanField(FIELD_NAME, i -> ((TestFieldItem<Boolean>) i).Value, isCompared);
    }

    @Test
    public void nameAndIsComparedReflectConstructorArgs()
    {
        assertEquals(FIELD_NAME, field(true).name());
        assertTrue(field(true).isCompared());
        assertFalse(field(false).isCompared());
    }

    @Test
    public void displayValueIsUpperCase()
    {
        BooleanField field = field(true);
        assertEquals("TRUE", field.displayValue(new TestFieldItem<>(true)));
        assertEquals("FALSE", field.displayValue(new TestFieldItem<>(false)));
    }

    @Test
    public void hasDiffIsFalseForEqualValues()
    {
        BooleanField field = field(true);
        assertFalse(field.hasDiff(new TestFieldItem<>(true), new TestFieldItem<>(true)));
        assertFalse(field.hasDiff(new TestFieldItem<>(false), new TestFieldItem<>(false)));
    }

    @Test
    public void hasDiffIsTrueForDifferentValues()
    {
        BooleanField field = field(true);
        assertTrue(field.hasDiff(new TestFieldItem<>(true), new TestFieldItem<>(false)));
        assertTrue(field.hasDiff(new TestFieldItem<>(false), new TestFieldItem<>(true)));
    }

    @Test
    public void determineDiffsFormatsValuesWhenDifferent()
    {
        BooleanField field = field(true);
        List<String> diffs = field.determineDiffs(new TestFieldItem<>(true), new TestFieldItem<>(false));
        assertEquals(List.of("BoolField(TRUE/FALSE)"), diffs);
    }

    @Test
    public void determineDiffsIsEmptyWhenEqual()
    {
        BooleanField field = field(true);
        assertTrue(field.determineDiffs(new TestFieldItem<>(true), new TestFieldItem<>(true)).isEmpty());
    }

    @Test
    public void withComparedReturnsNewFieldWithUpdatedIsCompared()
    {
        BooleanField field = field(true);
        Field updated = field.withCompared(false);
        assertFalse(updated.isCompared());
        assertEquals("TRUE", updated.displayValue(new TestFieldItem<>(true)));
    }

    @Test
    public void withThresholdsAreUnsupportedAndThrow()
    {
        BooleanField field = field(true);
        assertThrows(UnsupportedFieldOverrideException.class, () -> field.withAbsoluteThreshold(5.0));
        assertThrows(UnsupportedFieldOverrideException.class, () -> field.withPercentThreshold(0.2));
    }
}
