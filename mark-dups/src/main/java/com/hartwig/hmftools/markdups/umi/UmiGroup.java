package com.hartwig.hmftools.markdups.umi;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.markdups.common.Fragment;

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

    public static final int MAX_UMI_BASE_DIFF = 1;

    public static boolean exceedsUmiIdDiff(final String first, final String second)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > MAX_UMI_BASE_DIFF)
                    return true;
            }
        }

        return false;
    }

}
