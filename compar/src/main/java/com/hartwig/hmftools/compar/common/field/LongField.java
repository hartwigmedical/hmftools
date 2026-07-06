package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class LongField implements Field
{
    public final String name;
    public final Function<ComparableItem, Long> extractValue;
    public final boolean isCompared;
    public final Double absoluteThreshold;
    public final Double percentThreshold;

    public LongField(final String name, final Function<ComparableItem, Long> extractValue, final boolean isCompared,
            final Double absoluteThreshold, final Double percentThreshold)
    {
        this.name = name;
        this.extractValue = extractValue;
        this.isCompared = isCompared;
        this.absoluteThreshold = absoluteThreshold;
        this.percentThreshold = percentThreshold;
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
    public Double absoluteThreshold()
    {
        return absoluteThreshold;
    }

    @Override
    public Double percentThreshold()
    {
        return percentThreshold;
    }

    @Override
    public String displayValue(final ComparableItem item)
    {
        return extractValue.apply(item).toString();
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        long oldValue = extractValue.apply(oldItem);
        long newValue = extractValue.apply(newItem);

        if(oldValue == newValue)
        {
            return false;
        }

        long absDiff = abs(oldValue - newValue);
        double relDiff = (double) absDiff / max(abs(oldValue), abs(newValue));

        boolean satisfiesAbsDiff = absoluteThreshold == null || absDiff > absoluteThreshold;
        boolean satisfiesRelDiff = percentThreshold == null || relDiff > percentThreshold;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }
}
