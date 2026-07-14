package com.hartwig.hmftools.compar.common.field;

import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class DisplayOnlyField implements Field
{
    // Only for showing extra info. Cannot be compared.
    private final String name;
    private final Function<ComparableItem, String> extractValue;
    private final Function<ComparableItem, Boolean> hasValue;

    public final static String DISPLAY_TYPE = "display";

    public DisplayOnlyField(final String name, final Function<ComparableItem, String> extractValue,
            final Function<ComparableItem, Boolean> hasValue)
    {
        this.name = name;
        this.extractValue = extractValue;
        this.hasValue = hasValue;
    }

    @Override
    public String name()
    {
        return name;
    }

    @Override
    public boolean isCompared()
    {
        return false;
    }

    @Override
    public String type()
    {
        return DISPLAY_TYPE;
    }

    @Override
    public String displayValue(final ComparableItem item)
    {
        return hasValue(item) ? extractValue.apply(item) : "";
    }

    @Override
    public boolean hasValue(final ComparableItem item)
    {
        return hasValue.apply(item);
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        return false;
    }
}
