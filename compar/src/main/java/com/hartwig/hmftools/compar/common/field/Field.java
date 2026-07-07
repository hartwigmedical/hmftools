package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

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

    String type();

    String displayValue(ComparableItem item);

    boolean hasDiff(ComparableItem oldItem, ComparableItem newItem);

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
