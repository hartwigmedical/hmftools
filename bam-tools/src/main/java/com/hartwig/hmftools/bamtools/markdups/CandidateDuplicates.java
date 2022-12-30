package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class CandidateDuplicates
{
    // incomplete fragments (ie missing a mate read) with a matching fragment coordinate, and so candidates for being duplicates
    public final int Position;
    public final List<Fragment> Fragments;

    public CandidateDuplicates(int position, final Fragment fragment)
    {
        Fragments = Lists.newArrayList(fragment);
        Position = position;
    }

    public String toString()
    {
        return format("%d: fragments(%d)", Position, Fragments.size());
    }
}
