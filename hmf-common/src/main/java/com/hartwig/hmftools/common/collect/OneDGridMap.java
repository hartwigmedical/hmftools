package com.hartwig.hmftools.common.collect;

import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.function.BinaryOperator;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

public class OneDGridMap<T>
{
    private final SortedMap<Integer, T> mContainer;

    public OneDGridMap()
    {
        mContainer = Maps.newTreeMap();
    }

    @Nullable
    public T get(int x)
    {
        return mContainer.get(x);
    }

    public void put(int x, final T value)
    {
        mContainer.put(x, value);
    }

    public void merge(int x, final T value, final BinaryOperator<T> accumulator)
    {
        T currentValue = mContainer.get(x);
        if(currentValue == null)
        {
            mContainer.put(x, value);
            return;
        }

        T newValue = accumulator.apply(currentValue, value);
        mContainer.put(x, newValue);
    }

    public List<T> mergeValuesByDistance(int maxDistance, final BinaryOperator<T> accumulator)
    {
        List<T> mergedValues = Lists.newArrayList();
        T currentMergedValue = null;
        int lastX = 0;
        for(Map.Entry<Integer, T> xAndValue : mContainer.entrySet())
        {
            int x = xAndValue.getKey();
            T value = xAndValue.getValue();

            if(currentMergedValue == null)
            {
                currentMergedValue = value;
                lastX = x;
                continue;
            }

            if(x - lastX <= maxDistance)
            {
                currentMergedValue = accumulator.apply(currentMergedValue, value);
                lastX = x;
                continue;
            }

            mergedValues.add(currentMergedValue);
            currentMergedValue = value;
            lastX = x;
        }

        if(currentMergedValue != null)
            mergedValues.add(currentMergedValue);

        return mergedValues;
    }
}
