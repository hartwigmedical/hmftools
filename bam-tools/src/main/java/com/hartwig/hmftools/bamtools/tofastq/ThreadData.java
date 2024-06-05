package com.hartwig.hmftools.bamtools.tofastq;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

// a simple class to manage thread local data
public class ThreadData
{
    private final AtomicInteger mNextId = new AtomicInteger();
    private final ThreadLocal<FastqWriterCache> mThreadFastqWriterCache;
    private final ThreadLocal<PartitionReader> mThreadPartitionReader;

    private final List<FastqWriterCache> mFastqWriterCacheList = Collections.synchronizedList(new ArrayList<>());
    private final List<PartitionReader> mPartitionReaderList = Collections.synchronizedList(new ArrayList<>());

    public List<PartitionReader> getAllThreadPartitionReaders()
    {
        return mPartitionReaderList;
    }

    public ThreadData(final ToFastqConfig config, final RemoteReadHandler remoteReadHandler)
    {
        mThreadFastqWriterCache = ThreadLocal.withInitial(() -> {
            // we need to assign a unique id
            FastqWriterCache fastqWriterCache = new FastqWriterCache(config, String.format("t%d", mNextId.incrementAndGet()));
            mFastqWriterCacheList.add(fastqWriterCache);
            return fastqWriterCache;
        });

        // make a PartitionReader per thread
        mThreadPartitionReader = ThreadLocal.withInitial(() -> {
            PartitionReader partitionReader = new PartitionReader(config, getFastqWriterCache(), remoteReadHandler);
            mPartitionReaderList.add(partitionReader);
            return partitionReader;
        });
    }

    public PartitionReader getPartitionReader()
    {
        return mThreadPartitionReader.get();
    }

    public FastqWriterCache getFastqWriterCache()
    {
        return mThreadFastqWriterCache.get();
    }

    public List<FastqWriterCache> getAllThreadFastqWriterCaches()
    {
        return mFastqWriterCacheList;
    }

    public void closePartitionReaders()
    {
        for(PartitionReader partitionReader : mPartitionReaderList)
        {
            partitionReader.close();
        }
    }

    public void closeFastqWriters()
    {
        for(FastqWriterCache fastqWriterCache : mFastqWriterCacheList)
        {
            fastqWriterCache.close();
        }
    }
}
