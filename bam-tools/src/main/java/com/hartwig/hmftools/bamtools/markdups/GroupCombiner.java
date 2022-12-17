package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import htsjdk.samtools.SAMRecord;

public class GroupCombiner
{
    private final RecordWriter mRecordWriter;

    private final Map<String,Map<String,DuplicateStatus>> mReadGroupStatus; // keyed by chromosome-partition then readId
    private final Map<String,Map<String,List<SAMRecord>>> mSupplementaries; // keyed by chromosome-partition then readId
    private final Set<String> mProcessedPartitions;

    private static final String CHR_PARTITION_DELIM = "_";

    public GroupCombiner(final RecordWriter recordWriter)
    {
        mRecordWriter = recordWriter;
        mReadGroupStatus = Maps.newHashMap();
        mSupplementaries = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
    }

    public synchronized void partitionComplete(final String chrPartition)
    {
        mProcessedPartitions.add(chrPartition);

        // clear up any pending
        Map<String,List<SAMRecord>> suppMap = mSupplementaries.get(chrPartition);

        if(suppMap != null)
        {
            mSupplementaries.remove(chrPartition);

            for(List<SAMRecord> suppGroups : suppMap.values())
            {
                suppGroups.forEach(x -> mRecordWriter.writeRecord(x, DuplicateStatus.NONE));
            }
        }

        mReadGroupStatus.remove(chrPartition);

        if((mProcessedPartitions.size() % 100) == 0)
        {
            int cachedSuppCount = mSupplementaries.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();
            int cachedReadGroups = mReadGroupStatus.values().stream().mapToInt(x -> x.size()).sum();

            BM_LOGGER.info("group cache partitions processed({}) cached supplementaries({}) readGroups({})",
                    mProcessedPartitions.size(), cachedSuppCount, cachedReadGroups);
        }
    }

    public void handleRemaining()
    {
        for(Map<String,List<SAMRecord>> readMap : mSupplementaries.values())
        {
            for(List<SAMRecord> reads : readMap.values())
            {
                reads.forEach(x -> mRecordWriter.writeRecord(x, DuplicateStatus.NONE));
            }
        }
    }

    public synchronized void processReadGroup(final String readId, final List<String> remotePartitions, DuplicateStatus dupStatus)
    {
        for(String chrPartition : remotePartitions)
        {
            if(mProcessedPartitions.contains(chrPartition))
            {
                checkSupplementaries(chrPartition, readId, dupStatus);
            }
            else
            {
                // store status for when this partition is processed
                Map<String,DuplicateStatus> groupStatus = mReadGroupStatus.get(chrPartition);

                if(groupStatus == null)
                {
                    groupStatus = Maps.newHashMap();
                    mReadGroupStatus.put(chrPartition, groupStatus);
                }

                groupStatus.put(readId, dupStatus);
            }
        }
    }

    private void checkSupplementaries(final String chrPartition, final String readId, DuplicateStatus dupStatus)
    {
        Map<String,List<SAMRecord>> suppMap = mSupplementaries.get(chrPartition);

        if(suppMap == null)
            return;

        List<SAMRecord> supplementaries = suppMap.get(readId);

        if(supplementaries != null)
        {
            supplementaries.forEach(x -> mRecordWriter.writeRecord(x, dupStatus));
            suppMap.remove(readId);
        }
    }

    public synchronized void handleSupplementary(final SAMRecord read, final String chrPartition)
    {
        Map<String,List<SAMRecord>> suppMap = mSupplementaries.get(chrPartition);

        if(suppMap == null)
        {
            suppMap = Maps.newHashMap();
            mSupplementaries.put(chrPartition, suppMap);
        }

        List<SAMRecord> supplementaries = suppMap.get(read.getReadName());

        if(supplementaries == null)
        {
            supplementaries = Lists.newArrayList();
            suppMap.put(read.getReadName(), supplementaries);
        }

        supplementaries.add(read);
    }

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }
}
