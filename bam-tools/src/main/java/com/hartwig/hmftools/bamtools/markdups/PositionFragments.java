package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class PositionFragments
{
    public final int Position;
    public final List<Fragment> Fragments;

    public PositionFragments(final Fragment fragment)
    {
        this(fragment.initialPosition(), fragment);
    }

    public PositionFragments(int position, final Fragment fragment)
    {
        Fragments = Lists.newArrayList(fragment);
        Position = position;
    }

    public PositionFragments(int position, final List<Fragment> fragments)
    {
        Fragments = fragments;
        Position = position;
    }

    public String toString()
    {
        return format("%d: fragments(%d)", Position, Fragments.size());
    }
}
