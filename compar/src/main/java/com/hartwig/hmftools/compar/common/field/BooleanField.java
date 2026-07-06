package com.hartwig.hmftools.compar.common.field;

import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class BooleanField implements Field
{
    public final String name;
    public final Function<ComparableItem, Boolean> extractValue;
    public final boolean isCompared;

    public BooleanField(final String name, final Function<ComparableItem, Boolean> extractValue, final boolean isCompared)
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
        return extractValue.apply(item).toString().toUpperCase();
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        return extractValue.apply(oldItem) != extractValue.apply(newItem);
    }
}
