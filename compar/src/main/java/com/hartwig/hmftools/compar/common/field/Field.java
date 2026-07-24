package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.compar.ComparableItem;

public interface Field
{
    String DISPLAY_VALUE_DELIMITER = ";";

    String name();

    boolean isCompared();

    default Double absoluteThreshold()
    {
        return null;
    }

    default Double percentThreshold()
    {
        return null;
    }

    default Field withCompared(final boolean compared)
    {
        if(compared == isCompared())
        {
            return this;
        }

        throw new UnsupportedFieldOverrideException(this, "compared");
    }

    default Field withAbsoluteThreshold(final Double absoluteThreshold)
    {
        if(Objects.equals(absoluteThreshold, absoluteThreshold()))
        {
            return this;
        }

        throw new UnsupportedFieldOverrideException(this, "absolute threshold");
    }

    default Field withPercentThreshold(final Double percentThreshold)
    {
        if(Objects.equals(percentThreshold, percentThreshold()))
        {
            return this;
        }

        throw new UnsupportedFieldOverrideException(this, "percent threshold");
    }

    String type();

    String displayValue(ComparableItem item);

    boolean hasDiff(ComparableItem oldItem, ComparableItem newItem);

    default boolean hasValue(ComparableItem item)
    {
        return true;
    }

    default List<String> determineDiffs(ComparableItem oldItem, ComparableItem newItem)
    {
        if(hasDiff(oldItem, newItem))
        {
            return List.of(format("%s(%s/%s)", name(), displayValue(oldItem), displayValue(newItem)));
        }
        else
        {
            return Collections.emptyList();
        }
    }
}
