package com.hartwig.hmftools.markdups.consensus;

import static java.lang.String.format;

public class GroupIdGenerator
{
    private int mCurrentId;

    public GroupIdGenerator()
    {
        mCurrentId = 1;
    }

    private static final int MAX_GROUP_ID = 999999;

    public long currentId() { return mCurrentId; }
    public void reset() { mCurrentId = 1; }

    public synchronized String nextId()
    {
        if(mCurrentId >= MAX_GROUP_ID)
            reset();

        return format("%06d", mCurrentId++);
    }

}
