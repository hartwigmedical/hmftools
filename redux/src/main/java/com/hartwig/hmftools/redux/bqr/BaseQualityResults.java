package com.hartwig.hmftools.redux.bqr;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.redux.BqrKey;

public class BaseQualityResults
{
    private final Map<BqrKey,Integer> mCombinedQualityCounts;

    private PerformanceCounter mPerfCounter;
    private long mTotalReadsUsed;
    private int mTotalAltsFiltered;

    public BaseQualityResults()
    {
        mCombinedQualityCounts = Maps.newHashMap();
        mPerfCounter = null;
        mTotalAltsFiltered = 0;
        mTotalReadsUsed = 0;
    }

    public void clear()
    {
        mCombinedQualityCounts.clear();
    }

    public synchronized void addBaseQualityRegionCounter(final BqrRegionReader regionCounter)
    {
        for(BqrKeyCounter counter : regionCounter.getQualityCounts())
        {
            Integer count = mCombinedQualityCounts.get(counter.Key);
            mCombinedQualityCounts.put(counter.Key, count != null ? count + counter.count() : counter.count());
        }

        mTotalReadsUsed += regionCounter.totalReadsUsed();
        mTotalAltsFiltered += regionCounter.totalAltsFiltered();
    }

    public synchronized void addPerfCounter(final PerformanceCounter perfCounter)
    {
        if(mPerfCounter == null)
        {
            mPerfCounter = perfCounter;
        }
        else
        {
            mPerfCounter.merge(perfCounter);
        }
    }

    public Map<BqrKey,Integer> getCombinedQualityCounts() { return mCombinedQualityCounts; }

    public long totalReadsUsed() { return mTotalReadsUsed; }
    public int totalAltsFiltered() { return mTotalAltsFiltered; }

    public void logPerfStats()
    {
        if(mPerfCounter != null)
            mPerfCounter.logStats();
    }
}
