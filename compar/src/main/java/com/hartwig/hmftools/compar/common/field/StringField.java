package com.hartwig.hmftools.compar.common.field;

import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class StringField implements Field
{
    private final String name;
    private final Function<ComparableItem, String> extractValue;
    private final boolean isCompared;

    public StringField(final String name, final Function<ComparableItem, String> extractValue, final boolean isCompared)
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
    public Field withCompared(final boolean compared)
    {
        return new StringField(name, extractValue, compared);
    }

    @Override
    public String type()
    {
        return "string";
    }

    @Override
    public String displayValue(final ComparableItem item)
    {
        return extractValue.apply(item);
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        return !extractValue.apply(oldItem).equals(extractValue.apply(newItem));
    }
}
