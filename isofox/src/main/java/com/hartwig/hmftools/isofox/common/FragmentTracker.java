package com.hartwig.hmftools.isofox.common;

import java.util.Map;

import com.google.common.collect.Maps;

public class FragmentTracker
{
    private final Map<String,Object> mReadMap;

    public FragmentTracker()
    {
        mReadMap = Maps.newHashMap();
    }

    public boolean hasReadId(final String readId) { return mReadMap.containsKey(readId); }
    public boolean hasRead(final ReadRecord read) { return mReadMap.containsKey(read.Id); }

    public int readsCount() { return mReadMap.size(); }

    public boolean checkReadId(final String readId)
    {
        Integer count = (Integer) mReadMap.get(readId);

        if(count == null)
        {
            mReadMap.put(readId, 1);
            return false;
        }

        mReadMap.remove(readId);
        return true;
    }

    public ReadRecord checkRead(final ReadRecord read)
    {
        ReadRecord otherRead = (ReadRecord)mReadMap.get(read.Id);

        if(otherRead == null)
        {
            mReadMap.put(read.Id, read);
            return null;
        }

        mReadMap.remove(read.Id);
        return otherRead;
    }

    public Object checkRead(final String readId, final Object store)
    {
        Object otherStore = mReadMap.get(readId);

        if(otherStore == null)
        {
            mReadMap.put(readId, store);
            return null;
        }

        mReadMap.remove(readId);
        return otherStore;
    }

    public void clear()
    {
        mReadMap.clear();
    }
}
