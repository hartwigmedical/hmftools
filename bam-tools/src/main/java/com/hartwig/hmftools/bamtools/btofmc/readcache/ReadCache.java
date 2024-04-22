package com.hartwig.hmftools.bamtools.btofmc.readcache;

import static java.lang.Math.max;

import static com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig.BFQ_LOGGER;

import java.util.Map;
import java.util.function.BiConsumer;

import com.google.common.collect.Maps;
import com.google.inject.Inject;
import com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig;

import htsjdk.samtools.SAMRecord;

public class ReadCache implements ReadCacheInterface
{
    private final BamToFastqConfig mConfig;
    private final BiConsumer<SAMRecord, SAMRecord> mReadPairConsumer;
    private final Map<String, SAMRecord> mCache;

    private int mMaxCacheSize;

    @Inject
    public ReadCache(final BamToFastqConfig config, final BiConsumer<SAMRecord, SAMRecord> readPairConsumer)
    {
        mConfig = config;
        mReadPairConsumer = readPairConsumer;
        mCache = Maps.newHashMap();

        mMaxCacheSize = 0;
    }

    @Override
    public void flush()
    {
    }

    @Override
    public int size()
    {
        return mCache.size();
    }

    @Override
    public boolean isEmpty()
    {
        return mCache.isEmpty();
    }

    @Override
    public void add(final SAMRecord read)
    {
        try
        {
            SAMRecord readPair = mCache.get(read.getReadName());
            if(readPair == null)
            {
                mCache.put(read.getReadName(), read);
                return;
            }

            mCache.remove(read.getReadName());
            mReadPairConsumer.accept(read, readPair);
        }
        finally
        {
            mMaxCacheSize = max(mMaxCacheSize, mCache.size());
        }
    }

    @Override
    public void mergeStats(final ReadCacheInterface o)
    {
        // TODO: this is kind of ugly.
        if(!(o instanceof ReadCache))
        {
            throw new RuntimeException("Cannot merge stats of non ReadCache into ReadCache");
        }

        ReadCache other = (ReadCache) o;
        // TODO: Should be me max or add?
        // TODO: Test this and use mutation testing.
        mMaxCacheSize += other.mMaxCacheSize;
    }

    @Override
    public void logStats()
    {
        if(mConfig.PerfDebug)
        {
            BFQ_LOGGER.debug("ReadCache stats: maxCacheSize({})", mMaxCacheSize);
        }
    }
}
