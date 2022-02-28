package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.max;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class BaseQualityResults
{
    private final List<BaseQualityRegionCounter> mRegionCounters;
    private PerformanceCounter mPerfCounter;

    public BaseQualityResults()
    {
        mRegionCounters = Lists.newArrayList();
        mPerfCounter = null;
    }

    public void clear()
    {
        mRegionCounters.clear();
    }

    public synchronized void addBaseQualityRegionCounter(final BaseQualityRegionCounter bqCounter)
    {
        mRegionCounters.add(bqCounter);
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

    public List<BaseQualityRegionCounter> getRegionCounters() { return mRegionCounters; }

    public void logPerfStats()
    {
        if(mPerfCounter != null)
            mPerfCounter.logStats();
    }
}
