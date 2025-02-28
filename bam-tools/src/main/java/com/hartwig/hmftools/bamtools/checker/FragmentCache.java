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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMRecord;

public class FragmentCache
{
    private final Map<String,Fragment> mFragmentMap;
    private final SAMFileHeader mSamFileHeader;

    private final FragmentStats mStats;
    private long mLastLogCount;
    private static final long LOG_READ_DIFF = 1_000_000;

    public FragmentCache(final SAMFileHeader fileHeader)
    {
        mFragmentMap = Maps.newHashMap();
        mSamFileHeader = fileHeader;

        mLastLogCount = 0;
        mStats = new FragmentStats();
    }

    public FragmentStats stats() { return mStats; }

    public synchronized List<SAMRecord> handleIncompleteFragments(final Collection<Fragment> fragments)
    {
        List<SAMRecord> completeReads = Lists.newArrayList();

        for(Fragment fragment : fragments)
        {
            Fragment existingFragment = mFragmentMap.get(fragment.readId());

            if(existingFragment == null)
            {
                serialiseFragmentReads(fragment);
                mFragmentMap.put(fragment.readId(), fragment);

                ++mStats.TotalFragments;
                ++mStats.InterPartitionFragments;

                continue;
            }

            existingFragment.transfer(fragment);

            // keep the fragment's reads if one or both primaries are missing
            if(!existingFragment.hasPrimaryInfo())
            {
                serialiseFragmentReads(existingFragment);
                continue;
            }

            deserialiseFragmentReads(existingFragment);

            List<SAMRecord> fragCompleteReads = existingFragment.extractCompleteReads();

            completeReads.addAll(fragCompleteReads);

            // keep just its primary info if its is waiting on supplementaries only
            if(existingFragment.isComplete())
            {
                mFragmentMap.remove(existingFragment.readId());

                if(existingFragment.expectedSupplementaryCount() > 0)
                    ++mStats.FragmentsWithSupplementaries;
            }
            else
            {
                existingFragment.serialiseReads();
            }
        }

        long currentCacheCount = readCount();

        if(abs(currentCacheCount - mLastLogCount) > LOG_READ_DIFF)
        {
            BT_LOGGER.info("fragment cache reads({} -> {})", mLastLogCount, currentCacheCount);
            mLastLogCount = currentCacheCount;
        }

        mStats.MaxFragmentCount = max(mStats.MaxFragmentCount, mFragmentMap.size());

        return completeReads;
    }

    public List<SAMRecord> extractReads()
    {
        if(mFragmentMap.isEmpty())
            return Collections.emptyList();

        List<SAMRecord> reads = Lists.newArrayList();

        for(Fragment fragment : mFragmentMap.values())
        {
            deserialiseFragmentReads(fragment);
            reads.addAll(fragment.reads());

            if(fragment.expectedSupplementaryCount() > 0)
                ++mStats.FragmentsWithSupplementaries;
        }

        return reads;
    }

    private void serialiseFragmentReads(final Fragment fragment)
    {
        if(mSamFileHeader == null)
            return;

        fragment.serialiseReads();
    }

    private void deserialiseFragmentReads(final Fragment fragment)
    {
        if(mSamFileHeader == null)
            return;

        fragment.deserialiseReads(mSamFileHeader);
    }

    public void clear() { mFragmentMap.clear(); }

    public long readCount() { return mFragmentMap.values().stream().mapToInt(x -> x.readCount()).sum(); }

    @VisibleForTesting
    public Map<String,Fragment> fragmentMap() { return mFragmentMap; }

    public long fragmentCount() { return mFragmentMap.size(); }


    public String toString()
    {
        return format("fragment cache size(%d reads=%d)", mFragmentMap.size(), readCount());
    }
}
