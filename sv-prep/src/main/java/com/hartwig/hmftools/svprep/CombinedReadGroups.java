package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.BAM;
import static com.hartwig.hmftools.svprep.WriteType.READS;
import static com.hartwig.hmftools.svprep.reads.ExpectedRead.getExpectedReads;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.svprep.reads.ExpectedRead;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadGroupState;
import com.hartwig.hmftools.svprep.reads.ReadRecord;

import org.apache.commons.math3.analysis.function.Exp;

public class CombinedReadGroups
{
    private final SvConfig mConfig;
    private final Map<String,Map<String,ReadGroupState>> mIncompleteReadGroups; // keyed by chromosome-partition then readId

    private final Map<String,Map<String,List<ExpectedRead>>> mChrReadGroupReads; // keyed by chromosome then readId

    private final PerformanceCounter mPerfCounter;
    private int mMatchedGroupCount;
    private int mUnmatchedGroupCount;
    private final Set<String> mProcessedPartitions;
    private int mLastSnapshotCount;

    public CombinedReadGroups(final SvConfig config)
    {
        mConfig = config;
        mIncompleteReadGroups = Maps.newHashMap();
        mChrReadGroupReads = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mMatchedGroupCount = 0;
        mUnmatchedGroupCount = 0;
        mLastSnapshotCount = 0;
        mPerfCounter = new PerformanceCounter("ReadMerge");
    }

    // public void logPerfStats() { mPerfCounter.logStats(); }

    private static final String CHR_PARTITION_DELIM = "_";
    private static final int LOG_CACH_DIFF = 50000;

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }

    public static String chromosomeFromChromosomePartition(final String chrPartition)
    {
        return chrPartition.split(CHR_PARTITION_DELIM, 2)[0];
    }

    public synchronized void processSpanningReadGroups(final Map<String,ReadGroup> spanningGroups)
    {
        // look for reads which have already been found (and therefore written)
        for(Map.Entry<String,ReadGroup> entry : spanningGroups.entrySet())
        {
            ReadGroup readGroup = entry.getValue();

            for(ReadRecord read : readGroup.reads())
            {
                processNewRead(readGroup, read, null);
            }

            // also the expected reads
            List<ExpectedRead> expectedReads = ExpectedRead.getExpectedReads(readGroup);

            for(ExpectedRead read : expectedReads)
            {
                processNewRead(readGroup, null, read);
            }
        }
    }

    private void processNewRead(final ReadGroup readGroup, final ReadRecord actualRead, final ExpectedRead expectedRead)
    {
        String chromosome = actualRead != null ? actualRead.Chromosome : expectedRead.Chromosome;

        if(ignoreChromosome(chromosome))
            return;

        ExpectedRead read = expectedRead != null ? expectedRead : ExpectedRead.fromRead(actualRead);

        Map<String,List<ExpectedRead>> readGroupReads = mChrReadGroupReads.get(chromosome);

        if(readGroupReads == null)
        {
            readGroupReads = Maps.newHashMap();
            mChrReadGroupReads.put(chromosome, readGroupReads);
        }

        List<ExpectedRead> existingGroupReads = readGroupReads.get(readGroup.id());

        if(existingGroupReads == null)
        {
            existingGroupReads = Lists.newArrayList();
            readGroupReads.put(readGroup.id(), existingGroupReads);
            existingGroupReads.add(read);
        }
        else
        {
            ExpectedRead matchedRead = existingGroupReads.stream().filter(x -> x.matches(read)).findFirst().orElse(null);
            if(matchedRead != null)
            {
                if(matchedRead.Found && actualRead != null)
                    actualRead.setWritten();

                matchedRead.setExpectedMatchCount(readGroup.partitionCount() - 1);
                matchedRead.registerMatch();

                if(matchedRead.fullyMatched())
                {
                    existingGroupReads.remove(matchedRead);

                    if(existingGroupReads.isEmpty())
                        readGroupReads.remove(readGroup.id());
                }

                return;
            }

            // add this new read
            existingGroupReads.add(read);
        }
    }

    private boolean ignoreChromosome(final String chrPartition)
    {
        if(mConfig.SpecificChromosomes.isEmpty())
            return false;

        String chromosome = chromosomeFromChromosomePartition(chrPartition);
        return !mConfig.SpecificChromosomes.contains(chromosome);
    }

    /*
    public synchronized void addIncompleteReadGroup(final String chrPartition, final Map<String,ReadGroupState> newIncompleteGroups)
    {
        // still need to process this partition even if empty, so pick up any other partitions waiting on it
        mPerfCounter.start();

        mProcessedPartitions.add(chrPartition);

        int initTotalIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();
        int initChrIncomplete = newIncompleteGroups.size();
        int totalMergedGroups = 0;

        // new groups look up by their own chr-partition, but store themselves against their remote chr-parition(s)

        for(Map.Entry<String,ReadGroupState> entry : newIncompleteGroups.entrySet())
        {
            // String otherChrPartition = entry.getKey();
            ReadGroupState groupState = entry.getValue();

            for(String otherChrPartition : groupState.RemoteChromosomePartitions)
            {
                if(ignoreChromosome(otherChrPartition))
                    continue;

                if(mProcessedPartitions.contains(otherChrPartition))
                    groupState.ProcessedChromosomePartitions.add(otherChrPartition);

                // first looks to see if this chr-partition has been registered by other partition already
                Map<String,ReadGroupState> existingGroups = mIncompleteReadGroups.get(otherChrPartition);

                if(existingGroups != null)
                {
                    // if so it reconciles to two and checks the combined state
                    int previousMerged = mMatchedGroupCount;

                    List<ReadGroupState> unmatchedGroups = reconcileReadGroups(otherChrPartition, existingGroups, chrPartition, newIncompleteGroups);
                    unmatchedReadGroups.addAll(unmatchedGroups);

                    int mergedCount = mMatchedGroupCount - previousMerged;
                    totalMergedGroups += mergedCount;

                    SV_LOGGER.trace("combined chromosome partitions pair({} & {}) existing({}) new({}) merged({}) missed({})",
                            chrPartition, otherChrPartition, existingGroups.size(), newIncompleteGroups.size(), mergedCount,
                            unmatchedGroups.size());
                }
                else
                {
                    // remote partition has been processed but had no groups matching this new one
                    unmatchedReadGroups.addAll(newIncompleteGroups.values());

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
        }

        // now clear out any group waiting on this new partition
        for(Map<String,ReadGroupState> groupStateMap : mIncompleteReadGroups.values())
        {
            List<ReadGroupState> unmatchedGroups = groupStateMap.values().stream()
                    .filter(x -> x.RemoteChromosomePartitions.contains(chrPartition))
                    .collect(Collectors.toList());

            unmatchedGroups.forEach(x -> groupStateMap.remove(x.ReadId));
        }

        int newTotalIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        if(abs(mLastSnapshotCount - newTotalIncomplete) > LOG_CACH_DIFF)
        {
            mLastSnapshotCount = newTotalIncomplete;
            SV_LOGGER.info("completed partitions({}) incompleteCount({})", mProcessedPartitions.size(), newTotalIncomplete);
        }

        SV_LOGGER.debug("chromosomePartition({}) merged({}) partials chrPartition({}) total({} -> {}) unmergedGroups({})",
                chrPartition, totalMergedGroups, initChrIncomplete, initTotalIncomplete, newTotalIncomplete, unmatchedReadGroups.size());

        mPerfCounter.stop();

        mUnmatchedGroupCount += unmatchedReadGroups.size();

        return unmatchedReadGroups;
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
                ++mMatchedGroupCount;
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
        // unmatchedGroups.addAll(existingUnmatchedGroups);

        return unmatchedGroups;
    }
    */

    public void writeRemainingReadGroups(final ResultsWriter writer, final Set<WriteType> writeTypes)
    {
        if(!writeTypes.contains(BAM) && !writeTypes.contains(READS))
            return;

        int remainingCached = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        SV_LOGGER.info("spanning partition groups: matched({}) unmatched({}) cached({})",
                mMatchedGroupCount, mUnmatchedGroupCount, remainingCached);

        mPerfCounter.logStats();
    }
}
