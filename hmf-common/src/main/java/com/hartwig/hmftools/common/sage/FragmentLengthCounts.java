package com.hartwig.hmftools.common.sage;

import java.util.Map;

import com.google.common.collect.Maps;

public class FragmentLengthCounts
{
    private final Map<Integer,int[]> mLengthCounts;

    public static final int REF_COUNT = 0;
    public static final int ALT_COUNT = 1;
    public static final int MAX_LENGTH = 500;

    public FragmentLengthCounts()
    {
        mLengthCounts = Maps.newHashMap();
    }

    public Map<Integer,int[]> lengthCounts() { return mLengthCounts; }

    public void addLength(final int length, final boolean isAlt)
    {
        if(length <= 0 || length > MAX_LENGTH)
            return;

        int[] counts = mLengthCounts.get(length);

        if(counts == null)
        {
            counts = new int[] {0, 0};
            mLengthCounts.put(length, counts);
        }

        if(isAlt)
            ++counts[ALT_COUNT];
        else
            ++counts[REF_COUNT];
    }

    public void merge(final FragmentLengthCounts other)
    {
        for(Map.Entry<Integer,int[]> entry : other.lengthCounts().entrySet())
        {
            int[] counts = mLengthCounts.get(entry.getKey());

            if(counts == null)
            {
                counts = new int[] {0, 0};
                mLengthCounts.put(entry.getKey(), counts);
            }

            counts[REF_COUNT] += entry.getValue()[REF_COUNT];
            counts[ALT_COUNT] += entry.getValue()[ALT_COUNT];
        }
    }
}
