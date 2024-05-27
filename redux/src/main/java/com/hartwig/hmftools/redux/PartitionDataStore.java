package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.PartitionData;

public class PartitionDataStore
{
    private final Map<String, PartitionData> mPartitionDataMap;
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

    public List<PartitionData> partitions() {return mPartitionDataMap.values().stream().collect(Collectors.toList()); }

    public String toString() { return format("partitions(%d)", mPartitionDataMap.size()); }
}
