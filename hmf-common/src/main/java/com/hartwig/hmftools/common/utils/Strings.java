package com.hartwig.hmftools.common.utils;

import java.util.List;
import java.util.StringJoiner;

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

        final StringJoiner combinedStr = new StringJoiner(String.valueOf(delim));
        sourceList.forEach(x -> combinedStr.add(x));
        return sourceList.toString();
    }

    public static String reverseString(final String str)
    {
        return new StringBuilder(str).reverse().toString();
    }

}
