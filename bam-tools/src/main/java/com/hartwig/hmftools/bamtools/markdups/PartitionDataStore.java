package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class PartitionDataStore
{
    private final RecordWriter mRecordWriter;

    private final Map<String,PartitionData> mPartitionDataMap;
    private final Lock mLock;

    public PartitionDataStore(final RecordWriter recordWriter)
    {
        mRecordWriter = recordWriter;
        mPartitionDataMap = Maps.newHashMap();
        mLock = new ReentrantLock();
    }

    public PartitionData getOrCreatePartitionData(final String chrPartition)
    {
        PartitionData partitionCache = mPartitionDataMap.get(chrPartition);

        if(partitionCache == null)
        {
            partitionCache = new PartitionData(chrPartition);
            mPartitionDataMap.put(chrPartition, partitionCache);
        }

        return partitionCache;
    }

    public List<PartitionData> partitions() {return mPartitionDataMap.values().stream().collect(Collectors.toList()); }

    public synchronized void checkPartitions()
    {
        /*
        if(!mLock.tryLock())
            return;

        try
        {
            for(PartitionData partitionData : mPartitionDataMap.values())
            {
                List<Fragment> resolvedFragments = partitionData.extractResolvedFragments();
                mRecordWriter.writeFragments(resolvedFragments);
            }
        }
        finally
        {
            mLock.unlock();
        }
        */
    }

    public String toString() { return format("partitions(%d)", mPartitionDataMap.size()); }
}
