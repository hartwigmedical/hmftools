package com.hartwig.hmftools.bamtools.metrics;

import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class CombinedStats
{
    private final CoverageMetrics mCoverageMetrics;
    private final FragmentLengths mFragmentLengths;
    private final PerformanceCounter mPerfCounter;
    private long mTotalReads;
    private final FlagStats mFlagStats;

    public CombinedStats(int maxCoverage)
    {
        mCoverageMetrics = new CoverageMetrics(maxCoverage);
        mFragmentLengths = new FragmentLengths();
        mPerfCounter = new PerformanceCounter("Coverage");
        mTotalReads = 0;
        mFlagStats = new FlagStats();
    }

    public CoverageMetrics coverageMetrics() { return mCoverageMetrics; }
    public FragmentLengths fragmentLengths() { return mFragmentLengths; }

    public PerformanceCounter perfCounter() { return mPerfCounter; }
    public long totalReads() { return mTotalReads; }

    public FlagStats flagStats()
    {
        return mFlagStats;
    }

    public synchronized void addStats(
            final CoverageMetrics metrics, final FragmentLengths fragmentLengths, int totalReads, final PerformanceCounter perfCounter,
            final FlagStats flagStats)
    {
        mCoverageMetrics.merge(metrics);
        mFragmentLengths.merge(fragmentLengths);
        mTotalReads += totalReads;
        mPerfCounter.merge(perfCounter);
        mFlagStats.merge(flagStats);
    }
}
