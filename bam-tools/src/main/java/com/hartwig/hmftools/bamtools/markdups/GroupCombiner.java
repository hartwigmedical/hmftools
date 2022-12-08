package com.hartwig.hmftools.bamtools.markdups;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import htsjdk.samtools.SAMRecord;

public class GroupCombiner
{
    private final MarkDupsConfig mConfig;
    private final RecordWriter mRecordWriter;

    private final Map<String, Map<String,DuplicateGroup>> mDuplicateGroups; // keyed by chromosome-partition then readId
    private final Set<String> mProcessedPartitions;

    private static final String CHR_PARTITION_DELIM = "_";

    public GroupCombiner(final MarkDupsConfig config, final RecordWriter recordWriter)
    {
        mConfig = config;
        mRecordWriter = recordWriter;
        mDuplicateGroups = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
    }

    public synchronized void addDuplicateGroups(final List<DuplicateGroup> duplicateGroups)
    {


    }

    public synchronized void handleSupplementary(final SAMRecord record)
    {

    }


    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }

    public static String chrFromChrPartition(final String chrPartition)
    {
        return chrPartition.split(CHR_PARTITION_DELIM, 2)[0];
    }

    private String chrPartition(final String chromosome, int position) { return formChromosomePartition(chromosome, position, mConfig.PartitionSize); }

}
