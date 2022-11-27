package com.hartwig.hmftools.isofox.common;

public class FragmentTypeCounts
{
    private long[] mTypeCounts;

    public FragmentTypeCounts()
    {
        mTypeCounts = new long[FragmentType.values().length];
    }

    public long typeCount(final FragmentType type) { return mTypeCounts[type.ordinal()]; }
    public void addCount(final FragmentType type) { ++mTypeCounts[type.ordinal()]; }
    public void addCount(final FragmentType type, long count) { mTypeCounts[type.ordinal()] += count; }

    public void combine(final FragmentTypeCounts other)
    {
        for(FragmentType type : FragmentType.values())
        {
            mTypeCounts[type.ordinal()] += other.typeCount(type);
        }
    }

    public void clear()
    {
        for(int i = 0; i < mTypeCounts.length; ++i)
        {
            mTypeCounts[i] = 0;
        }
    }
}
