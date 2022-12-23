package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class PositionFragments
{
    public final int Position;
    public final List<Fragment> Fragments;

    public PositionFragments(final Fragment fragment, int position)
    {
        Fragments = Lists.newArrayList(fragment);
        Position = position;
    }

    public String toString()
    {
        return format("%d: reads(%d)", Position, Fragments.size());
    }
}
