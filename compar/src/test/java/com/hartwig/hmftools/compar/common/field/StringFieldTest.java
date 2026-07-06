package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class StringFieldTest
{
    private static final String FIELD_NAME = "StrField";

    private static StringField field(final boolean isCompared)
    {
        return new StringField(FIELD_NAME, i -> ((TestFieldItem<String>) i).Value, isCompared);
    }

    @Test
    public void nameAndIsComparedReflectConstructorArgs()
    {
        assertEquals(FIELD_NAME, field(true).name());
        assertTrue(field(true).isCompared());
        assertFalse(field(false).isCompared());
    }

    @Test
    public void displayValueReturnsRawString()
    {
        StringField field = field(true);
        assertEquals("hello", field.displayValue(new TestFieldItem<>("hello")));
        assertEquals("", field.displayValue(new TestFieldItem<>("")));
    }

    @Test
    public void hasDiffIsFalseForEqualStrings()
    {
        StringField field = field(true);
        assertFalse(field.hasDiff(new TestFieldItem<>("TEST"), new TestFieldItem<>("TEST")));
        assertFalse(field.hasDiff(new TestFieldItem<>(""), new TestFieldItem<>("")));
    }

    @Test
    public void hasDiffIsTrueForDifferentStrings()
    {
        StringField field = field(true);
        assertTrue(field.hasDiff(new TestFieldItem<>("TEST"), new TestFieldItem<>("test")));
        assertTrue(field.hasDiff(new TestFieldItem<>(""), new TestFieldItem<>("TEST")));
    }

    @Test
    public void determineDiffsFormatsValuesWhenDifferent()
    {
        StringField field = field(true);
        List<String> diffs = field.determineDiffs(new TestFieldItem<>("ABC"), new TestFieldItem<>("XYZ"));
        assertEquals(List.of("StrField(ABC/XYZ)"), diffs);
    }

    @Test
    public void determineDiffsIsEmptyWhenEqual()
    {
        StringField field = field(true);
        assertTrue(field.determineDiffs(new TestFieldItem<>("ABC"), new TestFieldItem<>("ABC")).isEmpty());
    }
}
