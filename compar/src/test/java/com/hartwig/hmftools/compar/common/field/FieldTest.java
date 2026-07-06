package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.compar.ComparableItem;

import org.junit.Test;

public class FieldTest
{
    private static class MinimalField implements Field
    {
        private final boolean mHasDiff;

        MinimalField(final boolean hasDiff)
        {
            mHasDiff = hasDiff;
        }

        @Override
        public String name()
        {
            return "MinimalField";
        }

        @Override
        public boolean isCompared()
        {
            return true;
        }

        @Override
        public String displayValue(final ComparableItem item)
        {
            return ((TestFieldItem<String>) item).Value;
        }

        @Override
        public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
        {
            return mHasDiff;
        }
    }

    @Test
    public void defaultThresholdsAreNull()
    {
        MinimalField field = new MinimalField(false);
        assertNull(field.absoluteThreshold());
        assertNull(field.percentThreshold());
    }

    @Test
    public void determineDiffsIsEmptyWhenNoDiff()
    {
        MinimalField field = new MinimalField(false);
        TestFieldItem<String> oldItem = new TestFieldItem<>("A");
        TestFieldItem<String> newItem = new TestFieldItem<>("B");

        assertTrue(field.determineDiffs(oldItem, newItem).isEmpty());
    }

    @Test
    public void determineDiffsFormatsFieldNameAndValuesWhenDiff()
    {
        MinimalField field = new MinimalField(true);
        TestFieldItem<String> oldItem = new TestFieldItem<>("A");
        TestFieldItem<String> newItem = new TestFieldItem<>("B");

        assertEquals(List.of("MinimalField(A/B)"), field.determineDiffs(oldItem, newItem));
    }
}
