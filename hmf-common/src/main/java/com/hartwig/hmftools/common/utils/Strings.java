package com.hartwig.hmftools.common.utils;

public final class Strings
{
    public static String appendStr(final String dest, final String source, char delim)
    {
        return dest.isEmpty() ? source : dest + delim + source;
    }

    public static String reverseString(final String str)
    {
        return new StringBuilder(str).reverse().toString();
    }

    public static String last(String str)
    {
        return str.substring(str.length() - 1);
    }
}
