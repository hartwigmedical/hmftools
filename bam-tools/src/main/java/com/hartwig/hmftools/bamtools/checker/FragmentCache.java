package com.hartwig.hmftools.bamtools.checker;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;

public class FragmentCache
{
    private final Map<String,Fragment> mFragmentMap;

    private long mLastLogCount;
    private static final long LOG_READ_DIFF = 1_000_000;

    public FragmentCache()
    {
        mFragmentMap = Maps.newHashMap();
        mLastLogCount = 0;
    }

    public synchronized List<SAMRecord> handleIncompleteFragments(final Collection<Fragment> fragments)
    {
        List<SAMRecord> completeReads = Lists.newArrayList();

        for(Fragment fragment : fragments)
        {
            Fragment existingFragment = mFragmentMap.get(fragment.readId());

            if(existingFragment == null)
            {
                mFragmentMap.put(fragment.readId(), fragment);
                continue;
            }

            fragment.reads().forEach(x -> existingFragment.addRead(x));

            List<SAMRecord> fragCompleteReads = existingFragment.extractCompleteReads();

            if(!fragCompleteReads.isEmpty())
            {
                completeReads.addAll(fragCompleteReads);

                if(existingFragment.isComplete())
                    mFragmentMap.remove(existingFragment.readId());
            }
        }

        long currentCacheCount = readCount();

        if(abs(currentCacheCount - mLastLogCount) > LOG_READ_DIFF)
        {
            BT_LOGGER.info("fragment cache reads({} -> {})", mLastLogCount, currentCacheCount);
            mLastLogCount = currentCacheCount;
        }

        return completeReads;
    }

    public List<SAMRecord> extractReads()
    {
        if(mFragmentMap.isEmpty())
            return Collections.emptyList();

        List<SAMRecord> reads = Lists.newArrayList();
        mFragmentMap.values().forEach(x -> reads.addAll(x.reads()));
        return reads;
    }

    public void clear() { mFragmentMap.clear(); }

    private long readCount() { return mFragmentMap.values().stream().mapToInt(x -> x.reads().size()).sum(); }

    public String toString()
    {
        return format("fragment cache size(%d reads=%d)", mFragmentMap.size(), readCount());
    }
}
