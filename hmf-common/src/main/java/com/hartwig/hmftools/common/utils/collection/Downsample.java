package com.hartwig.hmftools.common.utils.collection;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class Downsample
{
    @NotNull
    public static <T> List<T> downsample(int maxEntries, @NotNull List<T> input)
    {
        if(input.size() <= maxEntries)
        {
            return input;
        }

        long scale = Math.round(Math.ceil(1.0 * input.size() / maxEntries));
        final List<T> result = Lists.newArrayList();

        for(int i = 0; i < input.size(); i++)
        {
            if(i % scale == 0)
            {
                result.add(input.get(i));
            }
        }
        return result;
    }
}
