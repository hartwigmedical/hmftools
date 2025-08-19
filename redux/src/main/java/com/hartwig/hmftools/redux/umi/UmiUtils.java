package com.hartwig.hmftools.redux.umi;

import static com.hartwig.hmftools.redux.ReduxConstants.DEFAULT_MAX_UMI_BASE_DIFF;

import java.util.Comparator;

import com.hartwig.hmftools.redux.common.DuplicateGroup;

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
            if(first.readCount() < second.readCount())
                return 1;
            else if(first.readCount() > second.readCount())
                return -1;
            else
                return 0;
        }
    }

    public static String trimPolyGTail(final String umiId)
    {
        int tailLength;
        for(tailLength = 0; tailLength < umiId.length(); tailLength++)
        {
            if(umiId.charAt(umiId.length() - 1 - tailLength) != 'G')
                break;
        }

        return umiId.substring(0, umiId.length() - tailLength);
    }

    public static int polyGTailLength(final String umiId)
    {
        int tailLength;
        for(tailLength = 0; tailLength < umiId.length(); tailLength++)
        {
            if(umiId.charAt(umiId.length() - 1 - tailLength) != 'G')
                break;
        }

        return tailLength;
    }
}
