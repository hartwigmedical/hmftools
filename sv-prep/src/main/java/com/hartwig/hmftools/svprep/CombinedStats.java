package com.hartwig.hmftools.svprep;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.svprep.reads.PartitionStats;

public class CombinedStats
{
    public final PartitionStats ReadStats;

    public final List<PerformanceCounter> PerfCounters;

    public CombinedStats()
    {
        ReadStats = new PartitionStats();
        PerfCounters = Lists.newArrayList();
    }

    public synchronized void addPartitionStats(final PartitionStats stats)
    {
        ReadStats.add(stats);
    }

    public synchronized void addPerfCounters(final PerformanceCounter[] perfCounters)
    {
        if(PerfCounters.isEmpty())
        {
            Arrays.stream(perfCounters).forEach(x -> PerfCounters.add(x));
        }

        for(int i = 0; i < PerfCounters.size(); ++i)
        {
            PerfCounters.get(i).merge(perfCounters[i]);
        }
    }
}
