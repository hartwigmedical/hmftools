package com.hartwig.hmftools.compar.common.field;

import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class StringListField implements Field
{
    public final String name;
    public final Function<ComparableItem, List<String>> extractValue;
    public final boolean isCompared;

    public StringListField(final String name, final Function<ComparableItem, List<String>> extractValue, final boolean isCompared)
    {
        this.name = name;
        this.extractValue = extractValue;
        this.isCompared = isCompared;
    }

    @Override
    public String name()
    {
        return name;
    }

    @Override
    public boolean isCompared()
    {
        return isCompared;
    }

    @Override
    public String displayValue(final ComparableItem item)
    {
        return String.join(DISPLAY_VALUE_DELIMITER, extractValue.apply(item));
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        return !extractValue.apply(oldItem).equals(extractValue.apply(newItem));
    }
}
