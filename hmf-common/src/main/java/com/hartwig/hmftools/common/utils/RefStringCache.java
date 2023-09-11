package com.hartwig.hmftools.common.utils;

import static java.lang.String.format;

import java.util.HashMap;
import java.util.Map;

public class RefStringCache
{
    private final Map<String,String> mMap;
    private final int mMaxCacheSize;

    private static final int DEFAULT_MAX_CACHE = 100_000;

    public RefStringCache(int maxCacheSize)
    {
        mMap = new HashMap<>();
        mMaxCacheSize = maxCacheSize;
    }

    public RefStringCache()
    {
        this(DEFAULT_MAX_CACHE);
    }

    public String intern(final String s)
    {
        if(mMap.size() >= mMaxCacheSize)
            mMap.clear();

        String exist = mMap.putIfAbsent(s, s);

        return (exist == null) ? s : exist;
    }

    public int size() { return mMap.size(); }

    public String toString() { return format("cacheSize(%d)", mMap.size()); }
}
