package com.hartwig.hmftools.bamtools.btofmc.util;

import static java.lang.String.format;

import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.jetbrains.annotations.NotNull;
import org.pcollections.PSortedMap;
import org.pcollections.TreePMap;

// TODO NEXT: TEST
public class PerformanceCounterLock implements Lock
{
    private final String mName;
    private final Lock mLock;
    private final AtomicReference<PSortedMap<Long, PerformanceCounter>> mLockWaitPcByThreadRef;

    public PerformanceCounterLock(final String name)
    {
        mName = name;
        mLock = new ReentrantLock();
        mLockWaitPcByThreadRef = new AtomicReference<>(TreePMap.empty());
    }

    private PerformanceCounter getLockWaitPc()
    {
        long threadId = Thread.currentThread().getId();
        PerformanceCounter pc = mLockWaitPcByThreadRef.get().get(threadId);
        if(pc != null)
        {
            return pc;
        }

        return mLockWaitPcByThreadRef.updateAndGet(ref ->
        {
            if(ref.containsKey(threadId))
            {
                return ref;
            }

            PerformanceCounter newPc = new PerformanceCounter(format("lock(%s) threadId(%d) waiting", mName, threadId), false);
            return ref.plus(threadId, newPc);
        }).get(threadId);
    }

    public void mergePerformanceCounter(final PerformanceCounterLock other)
    {
        PSortedMap<Long, PerformanceCounter> otherLockWaitPcByThread = other.mLockWaitPcByThreadRef.get();
        mLockWaitPcByThreadRef.getAndUpdate(ref ->
        {
            for(Map.Entry<Long, PerformanceCounter> entry : otherLockWaitPcByThread.entrySet())
            {
                long threadId = entry.getKey();
                PerformanceCounter pc = entry.getValue();
                if(ref.containsKey(threadId))
                {
                    ref.get(threadId).merge(pc);
                    continue;
                }

                ref = ref.plus(threadId, pc);
            }

            return ref;
        });
    }

    public void logStats()
    {
        mLockWaitPcByThreadRef.get().values().forEach(PerformanceCounter::logStats);
    }

    @Override
    public void lock()
    {
        PerformanceCounter pc = getLockWaitPc();
        pc.start();
        mLock.lock();
        pc.stop();
    }

    @Override
    public void lockInterruptibly() throws InterruptedException
    {
        mLock.lockInterruptibly();
    }

    @Override
    public boolean tryLock()
    {
        return mLock.tryLock();
    }

    @Override
    public boolean tryLock(final long time, final TimeUnit unit) throws InterruptedException
    {
        return mLock.tryLock(time, unit);
    }

    @Override
    public void unlock()
    {
        mLock.unlock();
    }

    @NotNull
    @Override
    public Condition newCondition()
    {
        return mLock.newCondition();
    }
}
