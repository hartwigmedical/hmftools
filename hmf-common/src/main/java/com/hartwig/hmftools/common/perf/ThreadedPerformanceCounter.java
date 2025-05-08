package com.hartwig.hmftools.common.perf;

import static java.lang.String.format;

import java.util.concurrent.atomic.AtomicReference;

import org.pcollections.PSortedMap;
import org.pcollections.TreePMap;

public class ThreadedPerformanceCounter
{
    private final String mName;
    private final boolean mTrackTimes;
    private final AtomicReference<PSortedMap<Long, PerformanceCounter>> mThreadPc = new AtomicReference<>(TreePMap.empty());

    public ThreadedPerformanceCounter(final String name, boolean trackTimes)
    {
        mName = name;
        mTrackTimes = trackTimes;
    }

    public PerformanceCounter getPerformanceCounter()
    {
        final long threadId = Thread.currentThread().getId();
        return mThreadPc.updateAndGet(currentRef ->
        {
            if(currentRef.containsKey(threadId))
                return currentRef;

            return currentRef.plus(threadId, new PerformanceCounter(format("%s threadId(%d)", mName, threadId), mTrackTimes));
        }).get(threadId);
    }

    public void logStats()
    {
        mThreadPc.get().values().forEach(PerformanceCounter::logStats);
    }
}
