package com.hartwig.hmftools.common.collect;

import java.util.List;
import java.util.function.BinaryOperator;
import java.util.function.Supplier;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public class PartitionUtils
{
    public static <T> List<T> mergePartitions(final List<List<T>> partitions, final BinaryOperator<T> accumulator)
    {
        return mergePartitions(partitions, accumulator, null);
    }

    public static <T> List<T> mergePartitions(final List<List<T>> partitions, final BinaryOperator<T> accumulator,
            @Nullable final Supplier<T> identity)
    {
        List<T> mergedValues = Lists.newArrayList();
        for(List<T> partition : partitions)
        {
            if(partition.isEmpty() && identity == null)
                continue;

            if(partition.isEmpty())
            {
                mergedValues.add(identity.get());
                continue;
            }

            T firstValue = partition.get(0);
            T mergedValue = identity == null ? firstValue : accumulator.apply(identity.get(), firstValue);
            for(int i = 1; i < partition.size(); i++)
            {
                T value = partition.get(i);
                mergedValue = accumulator.apply(mergedValue, value);
            }

            mergedValues.add(mergedValue);
        }

        return mergedValues;
    }
}
