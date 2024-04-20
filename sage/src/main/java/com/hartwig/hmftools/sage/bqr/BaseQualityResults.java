package com.hartwig.hmftools.sage.bqr;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class BaseQualityResults
{
    private final Map<BqrKey,Integer> mCombinedQualityCounts;

    private PerformanceCounter mPerfCounter;

    public BaseQualityResults()
    {
        mCombinedQualityCounts = Maps.newHashMap();
        mPerfCounter = null;
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

    public void logPerfStats()
    {
        if(mPerfCounter != null)
            mPerfCounter.logStats();
    }
}
