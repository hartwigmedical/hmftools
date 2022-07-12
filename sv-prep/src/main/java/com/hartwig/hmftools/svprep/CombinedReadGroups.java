package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.BAM;
import static com.hartwig.hmftools.svprep.WriteType.READS;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.svprep.reads.ReadGroupState;

public class CombinedReadGroups
{
    private final SvConfig mConfig;
    private final Map<String,Map<String,ReadGroupState>> mIncompleteReadGroups; // keyed by chromosome then readId
    private final PerformanceCounter mPerfCounter;
    private int mMergedGroupCount;
    private int mUnmatchedGroupCount;
    private final Set<String> mProcessedPartitions;

    public CombinedReadGroups(final SvConfig config)
    {
        mConfig = config;
        mIncompleteReadGroups = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mMergedGroupCount = 0;
        mUnmatchedGroupCount = 0;
        mPerfCounter = new PerformanceCounter("ReadMerge");
    }

    // public void logPerfStats() { mPerfCounter.logStats(); }

    private static final String CHR_PARTITION_DELIM = "_";

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }

    public synchronized List<ReadGroupState> addIncompleteReadGroup(
            final String chrPartition, final Map<String,Map<String,ReadGroupState>> remoteChrIncompleteGroups)
    {
        List<ReadGroupState> unmatchedReadGroups = Lists.newArrayList();

        // still need to process this partition even if empty, so pick up any other partitions waiting on it
        //if(remoteChrIncompleteGroups.isEmpty())
        //    return unmatchedReadGroups;

        mPerfCounter.start();

        mProcessedPartitions.add(chrPartition);

        int initTotalIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();
        int initChrIncomplete = remoteChrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum();
        int totalMergedGroups = 0;

        // incomplete groups are looked up by the new chromosome partition to check for matches, and any which aren't found and
        // stored against the remote chromosome partition

        for(Map.Entry<String,Map<String,ReadGroupState>> entry : remoteChrIncompleteGroups.entrySet())
        {
            String otherChrPartition = entry.getKey();

            if(ignoreChromosome(otherChrPartition))
                continue;

            Map<String,ReadGroupState> newIncompleteGroups = entry.getValue();

            // first looks to see if this chr-partition has been registered by other partition already
            Map<String,ReadGroupState> existingGroups = mIncompleteReadGroups.get(otherChrPartition);

            if(existingGroups != null)
            {
                // if so it reconciles to two and checks the combined state
                int previousMerged = mMergedGroupCount;

                List<ReadGroupState> unmatchedGroups = reconcileReadGroups(chrPartition, existingGroups, newIncompleteGroups);
                unmatchedReadGroups.addAll(unmatchedGroups);

                int mergedCount = mMergedGroupCount - previousMerged;
                totalMergedGroups += mergedCount;

                SV_LOGGER.trace("combined chromosome partitions pair({} & {}) existing({}) new({}) merged({}) missed({})",
                        chrPartition, otherChrPartition, existingGroups.size(), newIncompleteGroups.size(), mergedCount,
                        unmatchedGroups.size());
            }
            else if(mProcessedPartitions.contains(otherChrPartition))
            {
                unmatchedReadGroups.addAll(newIncompleteGroups.values());
            }
            else
            {
                // otherwise the new chr-partition is cached against the one it come from
                Map<String,ReadGroupState> chrPartitionGroups = mIncompleteReadGroups.get(chrPartition);
                if(chrPartitionGroups == null)
                {
                    chrPartitionGroups = Maps.newHashMap();
                    mIncompleteReadGroups.put(chrPartition, chrPartitionGroups);
                }

                chrPartitionGroups.putAll(newIncompleteGroups);

                SV_LOGGER.trace("added chromosome partitions pair({} & {}) newGroups({})",
                        chrPartition, otherChrPartition, newIncompleteGroups.size());
            }
        }

        // now find any group waiting on this new partition
        for(Map<String,ReadGroupState> groupStateMap : mIncompleteReadGroups.values())
        {
            List<ReadGroupState> unmatchedGroups = groupStateMap.values().stream()
                    .filter(x -> x.RemoteChrPartition.equals(chrPartition))
                    .collect(Collectors.toList());

            unmatchedGroups.forEach(x -> groupStateMap.remove(x.ReadId));
            unmatchedReadGroups.addAll(unmatchedGroups);
        }

        int newTotalIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        SV_LOGGER.debug("chromosomePartition({}) merged({}) partials chrPartition({}) total({} -> {}) unmergedGroups({})",
                chrPartition, totalMergedGroups, initChrIncomplete, initTotalIncomplete, newTotalIncomplete, unmatchedReadGroups.size());

        mPerfCounter.stop();

        mUnmatchedGroupCount += unmatchedReadGroups.size();

        return unmatchedReadGroups;
    }

    private boolean ignoreChromosome(final String chrPartition)
    {
        if(mConfig.SpecificChromosomes.isEmpty())
            return false;

        String chromosome = chrPartition.split(CHR_PARTITION_DELIM, 2)[0];
        return !mConfig.SpecificChromosomes.contains(chromosome);
    }

    private List<ReadGroupState> reconcileReadGroups(
            final String newChrPartition, final Map<String,ReadGroupState> existingGroups, final Map<String,ReadGroupState> newGroups)
    {
        // all new groups and existing groups relate to the same partition, so check for matches on read Id
        // any unmatched groups where the pair of partitions match can then be purged since they will have no other chance to match
        List<ReadGroupState> unmatchedGroups = Lists.newArrayList();

        for(Map.Entry<String,ReadGroupState> entry : newGroups.entrySet())
        {
            // look for an existing incomplete group to add these reads to
            final String readId = entry.getKey();
            ReadGroupState existingGroupState = existingGroups.get(readId);

            if(existingGroupState != null)
            {
                // existingReadGroup.merge(srcGroupState); // no need for this
                ++mMergedGroupCount;
                existingGroups.remove(readId);
            }
            else
            {
                unmatchedGroups.add(entry.getValue());
            }
        }

        List<ReadGroupState> existingUnmatchedGroups = existingGroups.values().stream()
                .filter(x -> x.RemoteChrPartition.equals(newChrPartition))
                .collect(Collectors.toList());

        existingUnmatchedGroups.forEach(x -> existingGroups.remove(x.ReadId));
        unmatchedGroups.addAll(existingUnmatchedGroups);

        return unmatchedGroups;
    }

    public void writeRemainingReadGroups(final ResultsWriter writer, final Set<WriteType> writeTypes)
    {
        if(!writeTypes.contains(BAM) && !writeTypes.contains(READS))
            return;

        int remainingCached = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        SV_LOGGER.info("spanning partition groups: merged({}) unmatched({}) cached({})",
                mMergedGroupCount, mUnmatchedGroupCount, remainingCached);

        mPerfCounter.logStats();
    }
}
