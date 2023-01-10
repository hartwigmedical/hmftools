package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class UmiGroup
{
    public final String UmiId;
    public final List<Fragment> Fragments;

    public UmiGroup(final String umiId, final Fragment fragment)
    {
        UmiId = umiId;
        Fragments = Lists.newArrayList(fragment);
    }

    public int fragmentCount() { return Fragments.size(); }

    public String toString() { return format("id(%s) fragments(%d)", UmiId, Fragments.size()); }

}
