package com.hartwig.hmftools.isofox.common;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class FragmentTracker
{
    private final Map<String,Object> mReadMap;

    public FragmentTracker()
    {
        mReadMap = Maps.newHashMap();
    }

    public List<Object> getValues() { return mReadMap.values().stream().collect(Collectors.toList()); }

    public int readsCount() { return mReadMap.size(); }

    public Read checkRead(final Read read)
    {
        Read otherRead = (Read)mReadMap.remove(read.Id);

        if(otherRead != null)
            return otherRead;

        mReadMap.put(read.Id, read);
        return null;
    }

    public Object checkRead(final String readId, final Object store)
    {
        Object otherStore = mReadMap.remove(readId);

        if(otherStore != null)
            return otherStore;

        mReadMap.put(readId, store);
        return null;
    }

    public void clear()
    {
        mReadMap.clear();
    }
}
