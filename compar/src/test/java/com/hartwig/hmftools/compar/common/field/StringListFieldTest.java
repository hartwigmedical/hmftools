package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class StringListFieldTest
{
    private static final String FIELD_NAME = "StrListField";

    private static StringListField field(final boolean isCompared)
    {
        return new StringListField(FIELD_NAME, i -> ((TestFieldItem<List<String>>) i).Value, isCompared);
    }

    @Test
    public void nameAndIsComparedReflectConstructorArgs()
    {
        assertEquals(FIELD_NAME, field(true).name());
        assertTrue(field(true).isCompared());
        assertFalse(field(false).isCompared());
    }

    @Test
    public void displayValueJoinsWithDelimiter()
    {
        StringListField field = field(true);
        assertEquals("A;B;C", field.displayValue(new TestFieldItem<>(List.of("A", "B", "C"))));
        assertEquals("", field.displayValue(new TestFieldItem<>(List.of())));
    }

    @Test
    public void hasDiffIsFalseForEqualLists()
    {
        StringListField field = field(true);
        assertFalse(field.hasDiff(new TestFieldItem<>(List.of("A", "B")), new TestFieldItem<>(List.of("A", "B"))));
        assertFalse(field.hasDiff(new TestFieldItem<>(List.of()), new TestFieldItem<>(List.of())));
    }

    @Test
    public void hasDiffIsTrueForDifferentContent()
    {
        StringListField field = field(true);
        assertTrue(field.hasDiff(new TestFieldItem<>(List.of("A", "B")), new TestFieldItem<>(List.of("A"))));
        assertTrue(field.hasDiff(new TestFieldItem<>(List.of()), new TestFieldItem<>(List.of("A"))));
    }

    @Test
    public void hasDiffIsTrueForDifferentOrder()
    {
        StringListField field = field(true);
        assertTrue(field.hasDiff(new TestFieldItem<>(List.of("A", "B")), new TestFieldItem<>(List.of("B", "A"))));
    }

    @Test
    public void determineDiffsFormatsJoinedValuesWhenDifferent()
    {
        StringListField field = field(true);
        List<String> diffs = field.determineDiffs(new TestFieldItem<>(List.of("A", "B")), new TestFieldItem<>(List.of("A")));
        assertEquals(List.of("StrListField(A;B/A)"), diffs);
    }

    @Test
    public void determineDiffsIsEmptyWhenEqual()
    {
        StringListField field = field(true);
        assertTrue(field.determineDiffs(new TestFieldItem<>(List.of("A")), new TestFieldItem<>(List.of("A"))).isEmpty());
    }

    @Test
    public void withComparedReturnsNewFieldWithUpdatedIsCompared()
    {
        StringListField field = field(true);
        Field updated = field.withCompared(false);
        assertFalse(updated.isCompared());
        assertEquals("A;B", updated.displayValue(new TestFieldItem<>(List.of("A", "B"))));
    }

    @Test
    public void withThresholdsAreUnsupportedAndThrow()
    {
        StringListField field = field(true);
        assertThrows(UnsupportedFieldOverrideException.class, () -> field.withAbsoluteThreshold(5.0));
        assertThrows(UnsupportedFieldOverrideException.class, () -> field.withPercentThreshold(0.2));
    }
}
