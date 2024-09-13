package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.PartitionData.PARTITION_CACHE_CANDIDATE_DUP_GROUPS;
import static com.hartwig.hmftools.redux.common.PartitionData.PARTITION_CACHE_DUP_GROUPS;
import static com.hartwig.hmftools.redux.common.PartitionData.PARTITION_CACHE_DUP_GROUP_READS;
import static com.hartwig.hmftools.redux.common.PartitionData.PARTITION_CACHE_INCOMPLETE_FRAGS;
import static com.hartwig.hmftools.redux.common.PartitionData.PARTITION_CACHE_RESOLVED_STATUS;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.PartitionData;

public class PartitionDataStore
{
    private final Map<String,PartitionData> mPartitionDataMap;
    private final ReduxConfig mConfig;

    public PartitionDataStore(final ReduxConfig config)
    {
        mConfig = config;
        mPartitionDataMap = Maps.newHashMap();
    }

    public synchronized PartitionData getOrCreatePartitionData(final String chrPartition)
    {
        PartitionData partitionCache = mPartitionDataMap.get(chrPartition);

        if(partitionCache == null)
        {
            partitionCache = new PartitionData(chrPartition, mConfig);

            if(mConfig.PerfDebug && mConfig.Threads > 1)
                partitionCache.togglePerfChecks();

            mPartitionDataMap.put(chrPartition, partitionCache);
        }

        return partitionCache;
    }

    public synchronized void logTotalCacheSize()
    {
        int[] totalCachedCounts = new int[PARTITION_CACHE_DUP_GROUP_READS+1];

        for(PartitionData partitionData : mPartitionDataMap.values())
        {
            int[] cachedCounts = partitionData.collectCacheCounts();

            for(int i = 0; i < cachedCounts.length; ++i)
            {
                totalCachedCounts[i] += cachedCounts[i];
            }
        }

        RD_LOGGER.info("partition cache total counts: incomplete({}) candidateGroups({}) resolved({}) dupGroups({}) dupGroupReads({})",
                totalCachedCounts[PARTITION_CACHE_INCOMPLETE_FRAGS], totalCachedCounts[PARTITION_CACHE_CANDIDATE_DUP_GROUPS],
                totalCachedCounts[PARTITION_CACHE_RESOLVED_STATUS], totalCachedCounts[PARTITION_CACHE_DUP_GROUPS],
                totalCachedCounts[PARTITION_CACHE_DUP_GROUP_READS]);
    }

    public List<PartitionData> partitions() {return mPartitionDataMap.values().stream().collect(Collectors.toList()); }

    public String toString() { return format("partitions(%d)", mPartitionDataMap.size()); }
}
