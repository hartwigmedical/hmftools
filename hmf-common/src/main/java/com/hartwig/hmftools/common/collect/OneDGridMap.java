package com.hartwig.hmftools.common.collect;

import static com.hartwig.hmftools.common.collect.PartitionUtils.mergePartitions;

import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.function.BinaryOperator;
import java.util.function.Supplier;

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

    public List<List<T>> partitionValuesByDistance(int maxDistance)
    {
        List<List<T>> partitions = Lists.newArrayList();
        List<T> currentPartition = null;
        int lastX = 0;
        for(Map.Entry<Integer, T> xAndValue : mContainer.entrySet())
        {
            int x = xAndValue.getKey();
            T value = xAndValue.getValue();

            if(currentPartition == null)
            {
                currentPartition = Lists.newArrayList(value);
                lastX = x;
                continue;
            }

            if(x - lastX <= maxDistance)
            {
                currentPartition.add(value);
                lastX = x;
                continue;
            }

            partitions.add(currentPartition);
            currentPartition = Lists.newArrayList(value);
            lastX = x;
        }

        if(currentPartition != null)
            partitions.add(currentPartition);

        return partitions;
    }

    public List<T> mergeValuesByDistance(int maxDistance, final BinaryOperator<T> accumulator)
    {
        return mergeValuesByDistance(maxDistance, accumulator, null);
    }

    public List<T> mergeValuesByDistance(int maxDistance, final BinaryOperator<T> accumulator, @Nullable final Supplier<T> identity)
    {
        return mergePartitions(partitionValuesByDistance(maxDistance), accumulator, identity);
    }
}
