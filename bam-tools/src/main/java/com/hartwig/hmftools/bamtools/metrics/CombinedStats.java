package com.hartwig.hmftools.bamtools.metrics;

import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class CombinedStats
{
    private final Metrics mMetrics;
    private final PerformanceCounter mPerfCounter;
    private long mTotalReads;

    public CombinedStats(int maxCoverage)
    {
        mMetrics = new Metrics(maxCoverage);
        mPerfCounter = new PerformanceCounter("Coverage");
        mTotalReads = 0;
    }

    public Metrics metrics() { return mMetrics; }
    public PerformanceCounter perfCounter() { return mPerfCounter; }
    public long totalReads() { return mTotalReads; }

    public synchronized void addStats(final Metrics metrics, int totalReads, final PerformanceCounter perfCounter)
    {
        mMetrics.merge(metrics);
        mTotalReads += totalReads;
        mPerfCounter.merge(perfCounter);
    }
}
