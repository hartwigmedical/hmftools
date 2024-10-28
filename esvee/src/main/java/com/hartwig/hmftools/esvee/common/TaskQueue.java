package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.Queue;

import com.google.common.annotations.VisibleForTesting;

public class TaskQueue<E extends Object>
{
    private final Queue<E> mQueue;
    private final int mInitialCount;
    private int mRemainingCount;
    private final int mLogCount;
    private final String mItemType;

    public TaskQueue(final Queue<E> queue, final String itemType, final int logCount)
    {
        mQueue = queue;
        mInitialCount = queue.size();
        mRemainingCount = mInitialCount;
        mItemType = itemType;
        mLogCount = logCount;
    }

    public int initialCount() { return mInitialCount; }

    public E removeItem()
    {
        checkRemainingCount();

        return mQueue.remove();
    }

    private synchronized void checkRemainingCount()
    {
        --mRemainingCount;
        int processed = mInitialCount - mRemainingCount;

        if(mLogCount > 0 && (processed % mLogCount) == 0)
        {
            SV_LOGGER.debug("processed {} {}, remaining({})", processed, mItemType, mRemainingCount);
        }
    }

    @VisibleForTesting
    public TaskQueue(final Queue<E> queue)
    {
        this(queue, "", 0);
    }
}
