package com.hartwig.hmftools.esvee.util;

import java.util.Comparator;
import java.util.function.Function;

public enum NaturalSortComparator implements Comparator<String>
{
    INSTANCE;

    public static <T> Comparator<T> of(final Function<T, String> keyExtractor)
    {
        return Comparator.comparing(keyExtractor, INSTANCE);
    }

    @Override
    public int compare(final String left, final String right)
    {
        int leftIndex = 0, rightIndex = 0;
        while (true)
        {
            if(leftIndex == left.length() && rightIndex == right.length())
                return 0;
            if(leftIndex == left.length())
                return -1;
            if(rightIndex == right.length())
                return 1;

            final char l = left.charAt(leftIndex++);
            final char r = right.charAt(rightIndex++);
            if(Character.isDigit(l) && Character.isDigit(r))
            {
                final String leftNumberString = takeInteger(left, leftIndex - 1);
                final int leftNumber = Integer.parseInt(leftNumberString);
                final String rightNumberString = takeInteger(right, rightIndex - 1);
                final int rightNumber = Integer.parseInt(rightNumberString);

                final int compare = Integer.compare(leftNumber, rightNumber);
                if(compare != 0)
                    return compare;

                leftIndex += leftNumberString.length() - 1;
                rightIndex += rightNumberString.length() - 1;
            }
            else
            {
                final int compare = Character.compare(l, r);
                if(compare != 0)
                    return compare;
            }
        }
    }

    private boolean isDigit(final char c)
    {
        // Unlike Character.isDigit, we don't want to handle non-ascii.
        return c >= '0' && c <= '9';
    }

    private String takeInteger(final String s, final int startIndex)
    {
        final StringBuilder sb = new StringBuilder();
        for(int i = startIndex; i < s.length(); i++)
            if(isDigit(s.charAt(i)))
                sb.append(s.charAt(i));
            else
                break;
        return sb.toString();
    }
}
