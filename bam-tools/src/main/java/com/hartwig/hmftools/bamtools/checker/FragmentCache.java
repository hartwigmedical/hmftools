package com.hartwig.hmftools.bamtools.checker;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class FragmentCache
{
    private final Map<String,Fragment> mFragmentMap;
    private final SAMFileHeader mSamFileHeader;

    private final PerformanceCounter mPerfCounter;

    private final FragmentStats mStats;
    private long mLastLogCount;
    private static final long LOG_READ_DIFF = 1_000_000;

    public FragmentCache(final SAMFileHeader fileHeader)
    {
        mFragmentMap = Maps.newConcurrentMap();
        mSamFileHeader = fileHeader;

        mLastLogCount = 0;
        mStats = new FragmentStats();

        mPerfCounter = new PerformanceCounter("FragmentCache");
    }

    public FragmentStats stats() { return mStats; }

    public List<SAMRecord> handleIncompleteFragments(final Collection<Fragment> fragments)
    {
        // locks are required around the fragment map when adding or removing elements, and around the fragments themselves when adding
        // or extracting reads
        mPerfCounter.start();

        List<SAMRecord> completeReads = Lists.newArrayList();

        for(Fragment fragment : fragments)
        {
            Fragment existingFragment = mFragmentMap.get(fragment.readId());

            if(existingFragment == null)
            {
                fragment.serialiseReads();
                mFragmentMap.put(fragment.readId(), fragment);

                ++mStats.TotalFragments;
                ++mStats.InterPartitionFragments;

                continue;
            }

            boolean isComplete = existingFragment.mergeFragment(fragment, mSamFileHeader, completeReads);

            if(isComplete)
            {
                mFragmentMap.remove(existingFragment.readId());

                if(existingFragment.expectedSupplementaryCount() > 0)
                    ++mStats.FragmentsWithSupplementaries;
            }
        }

        long currentCacheCount = readCount();

        if(abs(currentCacheCount - mLastLogCount) > LOG_READ_DIFF)
        {
            BT_LOGGER.info("fragment cache reads({} -> {})", mLastLogCount, currentCacheCount);
            mLastLogCount = currentCacheCount;
        }

        mStats.MaxFragmentCount = max(mStats.MaxFragmentCount, mFragmentMap.size());

        mPerfCounter.stop();

        return completeReads;
    }

    public List<SAMRecord> extractReads()
    {
        if(mFragmentMap.isEmpty())
            return Collections.emptyList();

        List<SAMRecord> reads = Lists.newArrayList();

        for(Fragment fragment : mFragmentMap.values())
        {
            fragment.deserialiseReads(mSamFileHeader);
            reads.addAll(fragment.reads());

            if(fragment.expectedSupplementaryCount() > 0)
                ++mStats.FragmentsWithSupplementaries;
        }

        return reads;
    }

    public void clear() { mFragmentMap.clear(); }

    public long readCount() { return mFragmentMap.values().stream().mapToInt(x -> x.readCount()).sum(); }

    public void logFinalStats()
    {
        BT_LOGGER.debug("fragment cache stats: ", mStats.toString());
        mPerfCounter.logStats();
    }

    @VisibleForTesting
    public Map<String,Fragment> fragmentMap() { return mFragmentMap; }

    public long fragmentCount() { return mFragmentMap.size(); }

    public String toString()
    {
        return format("fragment cache size(%d reads=%d)", mFragmentMap.size(), readCount());
    }
}
