package com.hartwig.hmftools.bamtools.btofmc.readcache;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig.BFQ_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.ToIntFunction;

import com.google.common.annotations.VisibleForTesting;
import com.google.inject.Inject;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.SAMRecord;

// TODO: Retry testing EvictingReadCache.
public class EvictingReadCache implements ReadCacheInterface
{
    private static final int DEFAULT_INITIAL_CAPACITY = 1 << 22;

    private static final int GROWTH_FACTOR = 2;
    private static final float MAX_LOAD_FACTOR = 0.25f;

    private List<SAMRecord> mOverflowCache;
    private final BiConsumer<SAMRecord, SAMRecord> mPairConsumer;
    private final boolean mIsResizeable;
    private final ToIntFunction<SAMRecord> mSamRecordHashCode;

    private SAMRecord[] mBuckets;
    private int mBucketedReadCount;

    private int mMaxCacheSize;
    private float mMaxOverflowProportion;
    private long mAddCount;
    private long mEvictionCount;

    @Inject
    public EvictingReadCache(final BiConsumer<SAMRecord, SAMRecord> pairConsumer)
    {
        this(pairConsumer, DEFAULT_INITIAL_CAPACITY, true, record -> record.getReadName().hashCode());
    }

    public EvictingReadCache(final BiConsumer<SAMRecord, SAMRecord> pairConsumer, int initialCapacity)
    {
        this(pairConsumer, initialCapacity, true, record -> record.getReadName().hashCode());
    }

    @VisibleForTesting
    public EvictingReadCache(final BiConsumer<SAMRecord, SAMRecord> pairConsumer, int initialCapacity, boolean isResizeable)
    {
        this(pairConsumer, initialCapacity, isResizeable, record -> record.getReadName().hashCode());
    }

    @VisibleForTesting
    public EvictingReadCache(final BiConsumer<SAMRecord, SAMRecord> pairConsumer, int initialCapacity, boolean isResizeable,
            final ToIntFunction<SAMRecord> samRecordHashCode)
    {
        mOverflowCache = Lists.newArrayList();
        mPairConsumer = pairConsumer;
        mIsResizeable = isResizeable;
        mSamRecordHashCode = samRecordHashCode;

        mBuckets = new SAMRecord[initialCapacity];
        mBucketedReadCount = 0;

        mMaxCacheSize = 0;
        mMaxOverflowProportion = 0;
        mAddCount = 0;
        mEvictionCount = 0;
    }

    @Override
    public int size()
    {
        return mBucketedReadCount + mOverflowCache.size();
    }

    @Override
    public boolean isEmpty()
    {
        return size() == 0;
    }

    @Override
    public void add(final SAMRecord read)
    {
        add(read, false);
    }

    private void add(final SAMRecord read, boolean isResizing)
    {
        try
        {
            int bucketIndex = mSamRecordHashCode.applyAsInt(read) % mBuckets.length;
            while(bucketIndex < 0)
            {
                bucketIndex += mBuckets.length;
            }

            if(mBuckets[bucketIndex] == null)
            {
                mBuckets[bucketIndex] = read;
                ++mBucketedReadCount;
                if(!isResizing)
                {
                    ++mAddCount;
                }
                return;
            }

            SAMRecord cachedRead = mBuckets[bucketIndex];
            if(read.getReadName().equals(cachedRead.getReadName()))
            {
                mBuckets[bucketIndex] = null;
                --mBucketedReadCount;
                mPairConsumer.accept(read, cachedRead);
                return;
            }

            mBuckets[bucketIndex] = read;
            mOverflowCache.add(cachedRead);
            if(!isResizing)
            {
                ++mAddCount;
                ++mEvictionCount;
            }
        }
        finally
        {
            mMaxCacheSize = max(mMaxCacheSize, size());
            mMaxOverflowProportion = max(mMaxOverflowProportion, 1.0f * mOverflowCache.size() / size());
            if(!isResizing)
            {
                checkLoadFactor();
            }
        }
    }

    private void checkLoadFactor()
    {
        if(!mIsResizeable)
        {
            return;
        }

        float loadFactor = 1.0f * mBucketedReadCount / mBuckets.length;
        if(loadFactor <= MAX_LOAD_FACTOR)
        {
            return;
        }

        BFQ_LOGGER.debug("Resizing EvictingReadCache from {} buckets to {} buckets", mBuckets.length, mBuckets.length * GROWTH_FACTOR);

        SAMRecord[] oldBuckets = mBuckets;
        mBuckets = new SAMRecord[mBuckets.length * GROWTH_FACTOR];
        mBucketedReadCount = 0;

        for(SAMRecord read : oldBuckets)
        {
            if(read == null)
            {
                continue;
            }

            add(read, true);
        }
    }

    @Override
    public void flush()
    {
        if(mBucketedReadCount > 0)
        {
            for(int i = 0; i < mBuckets.length; ++i)
            {
                if(mBuckets[i] == null)
                {
                    continue;
                }

                mOverflowCache.add(mBuckets[i]);
                mBuckets[i] = null;
            }

            mBucketedReadCount = 0;
        }

        Collections.sort(mOverflowCache, Comparator.comparing(x -> x.getReadName()));
        List<SAMRecord> newOverflowCache = Lists.newArrayList();
        int i = 0;
        while(i < mOverflowCache.size())
        {
            SAMRecord current = mOverflowCache.get(i);
            if(i == mOverflowCache.size() - 1)
            {
                newOverflowCache.add(current);
                break;
            }

            SAMRecord next = mOverflowCache.get(i + 1);
            if(current.getReadName().equals(next.getReadName()))
            {
                mPairConsumer.accept(current, next);
                i += 2;
                continue;
            }

            newOverflowCache.add(current);
            ++i;
        }

        mOverflowCache = newOverflowCache;
    }

    @Override
    public void mergeStats(final ReadCacheInterface o)
    {
        // TODO: this is kind of ugly.
        if(!(o instanceof EvictingReadCache))
        {
            throw new RuntimeException("Cannot merge stats of non EvictingReadCache into EvictingReadCache");
        }

        EvictingReadCache other = (EvictingReadCache) o;

        // TODO: Do we want addition vs max
        // TODO: Test this and use mutation testing.
        mMaxCacheSize = max(mMaxCacheSize, other.mMaxCacheSize);
        mMaxOverflowProportion = max(mMaxOverflowProportion, other.mMaxOverflowProportion);
        mAddCount += other.mAddCount;
        mEvictionCount += other.mEvictionCount;
    }

    @Override
    public void logStats()
    {
        BFQ_LOGGER.info(format("EvictingReadCache stats: maxCacheSize(%d) maxOverflowProportion(%.2f%%) evictionRate(%.2f%%)", mMaxCacheSize,
                100.0f * mMaxOverflowProportion, 100.0f * mEvictionCount / mAddCount));
    }

    @VisibleForTesting
    public long addCount()
    {
        return mAddCount;
    }

    @VisibleForTesting
    public long evictionCount()
    {
        return mEvictionCount;
    }
}
