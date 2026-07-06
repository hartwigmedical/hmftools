package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class DoubleField implements Field
{
    public final String name;
    public final Function<ComparableItem, Double> extractValue;
    public final boolean isCompared;
    public final Double absoluteThreshold;
    public final Double percentThreshold;
    public final String formatString;

    public DoubleField(final String name, final Function<ComparableItem, Double> extractValue, final boolean isCompared,
            final Double absoluteThreshold, final Double percentThreshold, final String formatString)
    {
        this.name = name;
        this.extractValue = extractValue;
        this.isCompared = isCompared;
        this.absoluteThreshold = absoluteThreshold;
        this.percentThreshold = percentThreshold;
        this.formatString = formatString;
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
        return format(formatString,  extractValue.apply(item));
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        double oldValue = extractValue.apply(oldItem);
        double newValue = extractValue.apply(newItem);

        if(oldValue == newValue)
        {
            return false;
        }

        double absDiff = abs(oldValue - newValue);
        double relDiff = absDiff / max(abs(oldValue), abs(newValue));

        boolean satisfiesAbsDiff = absoluteThreshold == null || absDiff > absoluteThreshold;
        boolean satisfiesRelDiff = percentThreshold == null || relDiff > percentThreshold;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }
}
