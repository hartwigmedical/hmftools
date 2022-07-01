package com.hartwig.hmftools.svprep;

import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class CombinedStats
{
    public final PartitionStats ReadStats;

    public final PerformanceCounter PerfCounter;

    public CombinedStats()
    {
        ReadStats = new PartitionStats();
        PerfCounter = new PerformanceCounter("Total");
    }

    public synchronized void addPartitionStats(final PartitionStats stats)
    {
        ReadStats.add(stats);
    }

    public synchronized void addPerfCounters(final PerformanceCounter perfCounter)
    {
        PerfCounter.merge(perfCounter);
    }
}
