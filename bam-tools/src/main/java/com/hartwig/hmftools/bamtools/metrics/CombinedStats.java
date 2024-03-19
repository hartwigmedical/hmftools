package com.hartwig.hmftools.bamtools.metrics;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class CombinedStats
{
    private final CoverageMetrics mCoverageMetrics;
    private final FragmentLengths mFragmentLengths;
    private final ReadCounts mReadCounts;
    private final FlagStats mFlagStats;
    private final Map<Integer,Integer> mOffTargetOverlapCounts;

    private final PerformanceCounter mPerfCounter;

    public CombinedStats(int maxCoverage)
    {
        mCoverageMetrics = new CoverageMetrics(maxCoverage);
        mFragmentLengths = new FragmentLengths();
        mPerfCounter = new PerformanceCounter("Coverage");
        mReadCounts = new ReadCounts();
        mFlagStats = new FlagStats();
        mOffTargetOverlapCounts = Maps.newHashMap();
    }

    public CoverageMetrics coverageMetrics() { return mCoverageMetrics; }
    public FragmentLengths fragmentLengths() { return mFragmentLengths; }

    public PerformanceCounter perfCounter() { return mPerfCounter; }
    public ReadCounts readCounts() { return mReadCounts; }

    public FlagStats flagStats()
    {
        return mFlagStats;
    }

    public Map<Integer,Integer> offTargetOverlapCounts() { return mOffTargetOverlapCounts; }

    public synchronized void addStats(
            final CoverageMetrics metrics, final FragmentLengths fragmentLengths, final ReadCounts readCounts,
            final FlagStats flagStats, final Map<Integer,Integer> offTargetOverlapCounts, final PerformanceCounter perfCounter)
    {
        mCoverageMetrics.merge(metrics);
        mFragmentLengths.merge(fragmentLengths);
        mReadCounts.merge(readCounts);
        mPerfCounter.merge(perfCounter);
        mFlagStats.merge(flagStats);

        for(Map.Entry<Integer,Integer> entry : offTargetOverlapCounts.entrySet())
        {
            Integer frequency = mOffTargetOverlapCounts.get(entry.getKey());
            mOffTargetOverlapCounts.put(entry.getKey(), frequency != null ? frequency + entry.getValue() : entry.getValue());
        }
    }
}
