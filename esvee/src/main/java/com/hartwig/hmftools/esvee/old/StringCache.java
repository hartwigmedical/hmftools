package com.hartwig.hmftools.esvee.old;

import java.lang.ref.SoftReference;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public final class StringCache
{
    private static final String[] INTEGERS = new String[1000];
    private static SoftReference<Map<String, String>> mStrings = new SoftReference<>(new HashMap<>());

    // This is dumb but saves substantial memory.
    public static String of(final char c)
    {
        switch(c)
        {
            case 'A':
                return "A";
            case 'T':
                return "T";
            case 'C':
                return "C";
            case 'G':
                return "G";
        }
        return String.valueOf(c);
    }

    public static String of(final int i)
    {
        if(i < INTEGERS.length)
        {
            if(INTEGERS[i] == null)
            {
                INTEGERS[i] = String.valueOf(i);
            }

            return INTEGERS[i];
        }

        return String.valueOf(i);
    }

    // @Contract("null -> null; !null -> !null")
    public static String tryDedupe(final String s)
    {
        if(s == null)
        {
            return null;
        }

        Map<String, String> dedupeMap = mStrings.get();
        if(dedupeMap == null)
        {
            dedupeMap = new HashMap<>();
            mStrings = new SoftReference<>(dedupeMap);
        }

        if(dedupeMap.size() > 10_000)
        {
            dedupeMap.clear();
        }

        return Objects.requireNonNullElse(dedupeMap.putIfAbsent(s, s), s);
    }
}
