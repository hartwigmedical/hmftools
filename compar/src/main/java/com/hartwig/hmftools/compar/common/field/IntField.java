package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import java.util.function.Function;

import com.hartwig.hmftools.compar.ComparableItem;

public class IntField implements Field
{
    private final String name;
    private final Function<ComparableItem, Integer> extractValue;
    private final boolean isCompared;
    private final Double absoluteThreshold;
    private final Double percentThreshold;

    public IntField(final String name, final Function<ComparableItem, Integer> extractValue, final boolean isCompared,
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
    public Field withCompared(final boolean compared)
    {
        return new IntField(name, extractValue, compared, absoluteThreshold, percentThreshold);
    }

    @Override
    public Field withAbsoluteThreshold(final Double absoluteThreshold)
    {
        return new IntField(name, extractValue, isCompared, absoluteThreshold, percentThreshold);
    }

    @Override
    public Field withPercentThreshold(final Double percentThreshold)
    {
        return new IntField(name, extractValue, isCompared, absoluteThreshold, percentThreshold);
    }

    @Override
    public String type()
    {
        return "int";
    }

    @Override
    public String displayValue(final ComparableItem item)
    {
        return item.isValid() ? extractValue.apply(item).toString() : "";
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        int oldValue = extractValue.apply(oldItem);
        int newValue = extractValue.apply(newItem);

        if(oldValue == newValue)
        {
            return false;
        }

        int absDiff = abs(oldValue - newValue);
        double relDiff = (double) absDiff / max(abs(oldValue), abs(newValue));

        boolean satisfiesAbsDiff = absoluteThreshold == null || absDiff > absoluteThreshold;
        boolean satisfiesRelDiff = percentThreshold == null || relDiff > percentThreshold;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }
}
