package com.hartwig.hmftools.bamtools.markdups;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import htsjdk.samtools.SAMRecord;

public class GroupCombiner
{
    private final MarkDupsConfig mConfig;
    private final RecordWriter mRecordWriter;

    private final Map<String,Map<String,Boolean>> mReadGroupStatus; // keyed by chromosome-partition then readId
    private final Map<String,Map<String,List<SAMRecord>>> mSupplementaries; // keyed by chromosome-partition then readId
    private final Set<String> mProcessedPartitions;

    private static final String CHR_PARTITION_DELIM = "_";

    public GroupCombiner(final MarkDupsConfig config, final RecordWriter recordWriter)
    {
        mConfig = config;
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

        if(suppMap == null)
            return;

        mSupplementaries.remove(chrPartition);

        for(List<SAMRecord> suppGroups : suppMap.values())
        {
            suppGroups.forEach(x -> mRecordWriter.writeRecord(x, false));
        }
    }

    public void handleRemaining()
    {
        for(Map<String,List<SAMRecord>> readMap : mSupplementaries.values())
        {
            for(List<SAMRecord> reads : readMap.values())
            {
                reads.forEach(x -> mRecordWriter.writeRecord(x, false));
            }
        }
    }

    public synchronized void processReadGroup(final String readId, final List<String> remotePartitions, boolean isDuplicate)
    {
        for(String chrPartition : remotePartitions)
        {
            if(mProcessedPartitions.contains(chrPartition))
            {
                checkSupplementaries(chrPartition, readId, isDuplicate);
            }
            else
            {
                // store status for when this partition is processed
                Map<String, Boolean> groupStatus = mReadGroupStatus.get(chrPartition);

                if(groupStatus == null)
                    groupStatus = Maps.newHashMap();

                groupStatus.put(readId, isDuplicate);
            }
        }
    }

    private void checkSupplementaries(final String chrPartition, final String readId, boolean isDuplicate)
    {
        Map<String,List<SAMRecord>> suppMap = mSupplementaries.get(chrPartition);

        if(suppMap == null)
            return;

        List<SAMRecord> supplementaries = suppMap.get(readId);

        if(supplementaries != null)
        {
            supplementaries.forEach(x -> mRecordWriter.writeRecord(x, isDuplicate));
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

    /*
    public static String chrFromChrPartition(final String chrPartition)
    {
        return chrPartition.split(CHR_PARTITION_DELIM, 2)[0];
    }

    private String chrPartition(final String chromosome, int position) { return formChromosomePartition(chromosome, position, mConfig.PartitionSize); }
    */
}
