package com.hartwig.hmftools.common.perf;

import static java.lang.String.format;

import java.util.HashMap;
import java.util.Map;

// a class for reusing commonly allocated strings to save on memory
// an alternative to String.intern() since allows the cache to be cleared
public class StringCache
{
    private final Map<String,String> mMap;
    private final int mMaxCacheSize;

    private static final int DEFAULT_MAX_CACHE = 100_000;

    public StringCache(int maxCacheSize)
    {
        mMap = new HashMap<>();
        mMaxCacheSize = maxCacheSize;
    }

    public StringCache()
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
