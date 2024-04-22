package com.hartwig.hmftools.bamtools.btofmc.readcache;

import java.util.function.Supplier;

import com.google.inject.Inject;
import com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig;
import com.hartwig.hmftools.bamtools.btofmc.util.PerformanceCounterLock;

import htsjdk.samtools.SAMRecord;

// TODO NEXT: TEST
public class LockableReadCache implements ReadCacheInterface
{
    private final BamToFastqConfig mConfig;
    private final ReadCacheInterface mReadCache;
    private final PerformanceCounterLock mLock;

    @Inject
    public LockableReadCache(final BamToFastqConfig config, final Supplier<ReadCacheInterface> readCacheSupplier)
    {
        mConfig = config;
        mReadCache = readCacheSupplier.get();
        mLock = new PerformanceCounterLock("ConcurrentReadCache");
    }

    public void lock()
    {
        mLock.lock();
    }

    public void unlock()
    {
        mLock.unlock();
    }

    @Override
    public void flush()
    {
        mReadCache.flush();
    }

    @Override
    public int size()
    {
        return mReadCache.size();
    }

    @Override
    public boolean isEmpty()
    {
        return mReadCache.isEmpty();
    }

    @Override
    public void add(final SAMRecord read)
    {
        mReadCache.add(read);
    }

    @Override
    public void mergeStats(final ReadCacheInterface o)
    {
        // TODO: this is kind of ugly.
        if(!(o instanceof LockableReadCache))
        {
            throw new RuntimeException("Cannot merge stats of non LockableReadCache into LockableReadCache");
        }

        LockableReadCache other = (LockableReadCache) o;
        mReadCache.mergeStats(other.mReadCache);
        mLock.mergePerformanceCounter(other.mLock);
    }

    @Override
    public synchronized void logStats()
    {
        mReadCache.logStats();
        if(mConfig.PerfDebug)
        {
            mLock.logStats();
        }
    }
}
