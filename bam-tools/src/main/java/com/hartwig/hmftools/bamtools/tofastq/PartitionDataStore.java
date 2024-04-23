package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;

public class PartitionDataStore
{
    private final Map<String, PartitionData> mPartitionDataMap;
    private final FastqConfig mConfig;

    public PartitionDataStore(final FastqConfig config)
    {
        mConfig = config;
        mPartitionDataMap = Maps.newHashMap();
    }

    public synchronized PartitionData getOrCreatePartitionData(final String chrPartition)
    {
        PartitionData partitionCache = mPartitionDataMap.get(chrPartition);

        if(partitionCache == null)
        {
            partitionCache = new PartitionData(chrPartition);

            if(mConfig.PerfDebug && mConfig.Threads > 1)
                partitionCache.togglePerfChecks();

            mPartitionDataMap.put(chrPartition, partitionCache);
        }

        return partitionCache;
    }

    public void processUnmatchedReads(final FastqWriterCache writerCache)
    {
        // now merge into a single map
        Map<String,SAMRecord> readMap = Maps.newHashMap();
        int unmatchedReadsCount = 0;

        for(PartitionData partitionData : mPartitionDataMap.values())
        {
            List<SAMRecord> unmatchedReads = partitionData.unmatchedReads();
            unmatchedReadsCount += unmatchedReads.size();

            for(SAMRecord read : unmatchedReads)
            {
                SAMRecord mate = readMap.remove(read);

                if(mate != null)
                {
                    writerCache.processReadPair(read, mate, null);
                }
                else
                {
                    readMap.put(read.getReadName(), read);
                }
            }
        }

        if(unmatchedReadsCount > 0)
        {
            BT_LOGGER.info("partition cache unmatched reads({})", unmatchedReadsCount);
        }

        if(readMap.isEmpty())
            return;

        BT_LOGGER.warn("final unmatched read({})", readMap.size());

        for(SAMRecord read : readMap.values())
        {
            BT_LOGGER.warn("unmatched read: {}", readToString(read));
        }
    }

    public List<PartitionData> partitions() {return mPartitionDataMap.values().stream().collect(Collectors.toList()); }

    public String toString() { return format("partitions(%d)", mPartitionDataMap.size()); }
}
