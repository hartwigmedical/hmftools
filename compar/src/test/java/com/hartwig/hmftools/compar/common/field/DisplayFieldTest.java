package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class DisplayFieldTest
{
    private static final String FIELD_NAME = "DisplayField";

    private static DisplayField field()
    {
        return new DisplayField(FIELD_NAME,
                i -> ((TestFieldItem<String>) i).Value,
                i -> ((TestFieldItem<String>) i).Value != null);
    }

    @Test
    public void nameAndTypeReflectConstructorArgsAndIsNeverCompared()
    {
        DisplayField field = field();
        assertEquals(FIELD_NAME, field.name());
        assertEquals(DisplayField.DISPLAY_TYPE, field.type());
        assertFalse(field.isCompared());
    }

    @Test
    public void hasValueDelegatesToConstructorFunction()
    {
        DisplayField field = field();
        assertTrue(field.hasValue(new TestFieldItem<>("value")));
        assertFalse(field.hasValue(new TestFieldItem<>(null)));
    }

    @Test
    public void displayValueReturnsExtractedValueWhenPresent()
    {
        DisplayField field = field();
        assertEquals("value", field.displayValue(new TestFieldItem<>("value")));
    }

    @Test
    public void displayValueIsEmptyWhenNoValue()
    {
        DisplayField field = field();
        assertEquals("", field.displayValue(new TestFieldItem<>(null)));
    }

    @Test
    public void hasDiffIsAlwaysFalse()
    {
        DisplayField field = field();
        assertFalse(field.hasDiff(new TestFieldItem<>("A"), new TestFieldItem<>("A")));
        assertFalse(field.hasDiff(new TestFieldItem<>("A"), new TestFieldItem<>("B")));
        assertFalse(field.hasDiff(new TestFieldItem<>(null), new TestFieldItem<>("B")));
    }

    @Test
    public void determineDiffsIsAlwaysEmpty()
    {
        DisplayField field = field();
        assertTrue(field.determineDiffs(new TestFieldItem<>("A"), new TestFieldItem<>("B")).isEmpty());
    }

    @Test
    public void withOverridesAreAllUnsupportedAndThrow()
    {
        DisplayField field = field();
        assertThrows(UnsupportedFieldOverrideException.class, () -> field.withCompared(true));
        assertThrows(UnsupportedFieldOverrideException.class, () -> field.withAbsoluteThreshold(5.0));
        assertThrows(UnsupportedFieldOverrideException.class, () -> field.withPercentThreshold(0.2));
    }
}
