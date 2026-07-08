package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CommonUtils;

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
        CommonUtils.warnUnsupportedFieldOverride(this, "compared");
        return this;
    }

    default Field withAbsoluteThreshold(final Double absoluteThreshold)
    {
        CommonUtils.warnUnsupportedFieldOverride(this, "absolute threshold");
        return this;
    }

    default Field withPercentThreshold(final Double percentThreshold)
    {
        CommonUtils.warnUnsupportedFieldOverride(this, "percent threshold");
        return this;
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
