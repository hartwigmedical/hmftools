package com.hartwig.hmftools.common.utils;

import java.util.List;

public final class Strings
{
    public static String appendStr(final String dest, final String source, char delim)
    {
        return dest.isEmpty() ? source : dest + delim + source;
    }

    public static String appendStrList(final List<String> sourceList, char delim)
    {
        if(sourceList.isEmpty())
            return "";

        final StringBuilder combinedStr = new StringBuilder(sourceList.get(0));

        for(int i = 1; i < sourceList.size(); ++i)
        {
            combinedStr.append(delim + sourceList.get(i));
        }

        return combinedStr.toString();
    }

    public static String reverseString(final String str)
    {
        String reverse = "";

        for(int i = str.length() - 1; i >= 0; --i)
        {
            reverse += str.charAt(i);
        }

        return reverse;
    }

}
