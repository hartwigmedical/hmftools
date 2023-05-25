package com.hartwig.hmftools.markdups.umi;

import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_MAX_UMI_BASE_DIFF;
import static com.hartwig.hmftools.markdups.common.FragmentCoordinates.FRAGMENT_REVERSED_ID;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.markdups.common.DuplicateGroup;
import com.hartwig.hmftools.markdups.common.Fragment;

public final class UmiUtils
{
    public static boolean exceedsUmiIdDiff(final String first, final String second)
    {
        return exceedsUmiIdDiff(first, second, DEFAULT_MAX_UMI_BASE_DIFF);
    }

    public static boolean exceedsUmiIdDiff(final String first, final String second, int permittedDiff)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > permittedDiff)
                    return true;
            }
        }

        return false;
    }

    public static int calcUmiIdDiff(final String first, final String second)
    {
        if(first.length() != second.length())
            return 0;

        int diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;
            }
        }

        return diffs;
    }

    public static class SizeComparator implements Comparator<DuplicateGroup>
    {
        public int compare(final DuplicateGroup first, final DuplicateGroup second)
        {
            if(first.fragments().size() < second.fragments().size())
                return 1;
            else if(first.fragments().size() > second.fragments().size())
                return -1;
            else
                return 0;
        }
    }

}
