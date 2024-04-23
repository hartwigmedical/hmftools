package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.NANO_IN_MILLISECOND;

import java.util.List;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;

public class PartitionData
{
    private final String mChrPartition;

    private final Map<String,SAMRecord> mUnmatchedPairs;

    // any update to the maps is done under a lock
    private Lock mLock;
    private long mLastCacheCount;
    private long mLockAcquireTime;
    private boolean mPerfChecks;

    private static final int LOG_CACHE_COUNT = 250000;

    public PartitionData(final String chrPartition)
    {
        mChrPartition = chrPartition;
        mUnmatchedPairs = Maps.newHashMap();

        mLock = new ReentrantLock();
        mLockAcquireTime = 0;
        mPerfChecks = false;
    }

    public List<SAMRecord> unmatchedReads() { return mUnmatchedPairs.values().stream().collect(Collectors.toList()); }

    private static final String CHR_PARTITION_DELIM = "_";

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosomeIndicator(chromosome) + partition;
    }

    public static String chromosomeIndicator(final String chromosome)
    {
        return chromosome + CHR_PARTITION_DELIM;
    }

    public List<ReadPair> processUnpairedReads(final List<SAMRecord> unmatchedReads)
    {
        try
        {
            acquireLock();

            List<ReadPair> readPairs = Lists.newArrayList();

            for(SAMRecord read : unmatchedReads)
            {
                SAMRecord mate = mUnmatchedPairs.remove(read.getReadName());

                if(mate != null)
                {
                    readPairs.add(new ReadPair(read, mate));
                }
                else
                {
                    mUnmatchedPairs.put(read.getReadName(), read);
                }
            }

            return readPairs;
        }
        finally
        {
            checkCachedCounts();
            mLock.unlock();
        }
    }

    private void checkCachedCounts()
    {
        long cacheCount = mUnmatchedPairs.size();

        if(abs(mLastCacheCount - cacheCount) < LOG_CACHE_COUNT)
            return;

        mLastCacheCount = cacheCount;

        BT_LOGGER.debug("partition({}) cached unmatchedReads({})", mChrPartition, mUnmatchedPairs.size());
    }

    private void acquireLock()
    {
        if(!mPerfChecks)
        {
            mLock.lock();
            return;
        }

        long startTime = System.nanoTime();
        mLock.lock();
        mLockAcquireTime += System.nanoTime() - startTime;
    }

    public void togglePerfChecks() { mPerfChecks = true; }
    public double totalLockTimeMs() { return mLockAcquireTime / NANO_IN_MILLISECOND; }

    public String toString()
    {
        return format("%s: status(%d) unmatchedReads(%d)", mChrPartition, mUnmatchedPairs.size());
    }

    @VisibleForTesting
    public Map<String,SAMRecord> unmatchedReadsMap() { return mUnmatchedPairs; }

    @VisibleForTesting
    public void clearState()
    {
        mUnmatchedPairs.clear();
    }
}
