package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;

public class SpanningReadCache
{
    private final SvConfig mConfig;
    private final int mPartitionSize;

    // a cache of read groups found from each chr-partition to help retrieve complete fragments (ie expected reads within a group)
    // when a new chr-partition completes, the following steps are done:
    // - record that the chr-partition has been processed (in mProcessedPartitions)
    // - search for any expected reads (ie mates from other read groups) in mChrPartitionReadGroupReads
    // - if a read is matched, make note of it has already been written (to BAM) or record that it has been for future matches
    // - remove any group once all reads have been matched/found
    // - for newly registered expected reads (ie in future chr-partitions), store the in mChrPartitionReadGroupReads
    // - also make a link of other (ie a 3rd) future chr-partitions back to the new chr-partition to speed up checking
    // - purge any unmatched groups from this chr-partition if they weren't matched

    /*
    - JUNCTIONs
        - if has no unprocessed remote partitions -> DONE
        - otherwise register readId only against remote partitions to be picked up as an expected readId
        - pick up any SUPPORT reads
            - add to read group and they will then be written
            - remove these from CombinedCache

    - SUPPORTs
        - use the expected ReadIds to pick these up to avoid passing to the cache at all
        - logically (due to synchronisation) will be either EXPECTED or have remote unprocessed partitions
        - if no unprocessed partitions then drop immediately -> DONE
        - otherwise cache as now (ie source partition and readId)

     */

    private final Map<String,Map<String,List<ReadRecord>>> mChrPartitionReadGroupReads; // keyed by chromosome-partition then readId
    private final Map<String,Map<String, CachedReadGroup>> mCandidatePartitionGroups; // keyed by chromosome-partition then readId

    // a map of remote chr-partition to read partition to readIds
    private final Map<String,Map<String,Set<String>>> mExpectedChrPartitionReadIds;

    private final Map<String,Set<String>> mJunctionPartitionReadIds;

    private final PerformanceCounter mPerfCounter;
    private final Set<String> mProcessedPartitions;

    private int mLastSnapshotCount;
    private int mMatchedCandidates;
    private int mPurgedCandidates;

    public SpanningReadCache(final SvConfig config)
    {
        mConfig = config;
        mPartitionSize = config.PartitionSize;
        mChrPartitionReadGroupReads = Maps.newHashMap();
        mCandidatePartitionGroups = Maps.newHashMap();
        mExpectedChrPartitionReadIds = Maps.newHashMap();
        mJunctionPartitionReadIds = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mLastSnapshotCount = 0;
        mMatchedCandidates = 0;
        mPurgedCandidates = 0;
        mPerfCounter = new PerformanceCounter("SpanningReads");
    }

    private static final String CHR_PARTITION_DELIM = "_";
    private static final int LOG_CACH_DIFF = 50000;

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }

    private static String chrFromChrPartition(final String chrPartition)
    {
        return chrPartition.split(CHR_PARTITION_DELIM, 2)[0];
    }

    private String chrPartition(final String chromosome, int position) { return formChromosomePartition(chromosome, position, mPartitionSize); }

    public synchronized Set<String> getExpectedReadIds(final ChrBaseRegion partitionRegion)
    {
        String chrPartition = chrPartition(partitionRegion.Chromosome, partitionRegion.start());

        Set<String> expectedReadIds = mJunctionPartitionReadIds.get(chrPartition);

        return expectedReadIds != null ? expectedReadIds : Sets.newHashSet();
    }

    public synchronized void processSpanningReadGroups(final ChrBaseRegion partitionRegion, final Map<String,ReadGroup> spanningGroups)
    {
        mPerfCounter.start();

        String sourceChrPartition = chrPartition(partitionRegion.Chromosome, partitionRegion.start());
        mProcessedPartitions.add(sourceChrPartition);

        // look for reads which have already been found (and therefore written)
        for(Map.Entry<String,ReadGroup> entry : spanningGroups.entrySet())
        {
            ReadGroup readGroup = entry.getValue();

            List<String> unprocessedPartitions = readGroup.remotePartitions().stream()
                    .filter(x -> !mProcessedPartitions.contains(x))
                    .filter(x -> mConfig.SpecificChromosomes.isEmpty() || !ignoreChromosome(chrFromChrPartition(x)))
                    .collect(Collectors.toList());

            if(readGroup.conditionalOnRemoteReads())
            {
                for(ReadRecord read : readGroup.reads())
                {
                    processCandidateRead(unprocessedPartitions, readGroup, read);
                }
            }
            else
            {
                processJunctionRead(unprocessedPartitions, readGroup, sourceChrPartition);
            }
        }

        // purge any expected reads which were already found
        purgePartition(sourceChrPartition);

        logCacheCount(false);

        mPerfCounter.stop();
    }

    private void processJunctionRead(
            final List<String> unprocessedPartitions, final ReadGroup readGroup, final String sourceChrPartition)
    {
        final ReadRecord read = readGroup.reads().get(0);

        if(ignoreChromosome(read.Chromosome))
            return;

        Map<String,CachedReadGroup> cachedReadGroups = mCandidatePartitionGroups.get(sourceChrPartition);

        if(cachedReadGroups != null)
        {
            CachedReadGroup cachedReadGroup = cachedReadGroups.get(readGroup.id());

            if(cachedReadGroup != null)
            {
                mMatchedCandidates += cachedReadGroup.Reads.size();
                cachedReadGroup.Reads.forEach(x -> readGroup.addRead(x));
                cachedReadGroup.Reads.clear();

                cachedReadGroups.remove(readGroup.id());

                // also purge from other remote partitions
                for(String otherRemotePartitions : cachedReadGroup.Partitions)
                {
                    if(!otherRemotePartitions.equals(sourceChrPartition))
                    {
                        Map<String, CachedReadGroup> otherReadGroups = mCandidatePartitionGroups.get(otherRemotePartitions);

                        if(otherReadGroups == null)
                            continue;

                        otherReadGroups.remove(readGroup.id());
                    }
                }
            }
        }

        // only store this junction readId if other remote unprocessed partitions within this same read group are expected
        if(unprocessedPartitions.isEmpty())
            return;

        // store the junction group's readId against each unprocessed partition, to use to capture expected reads
        for(String unprocessedPartition : unprocessedPartitions)
        {
            Set<String> readIds = mJunctionPartitionReadIds.get(unprocessedPartition);

            if(readIds == null)
            {
                readIds = Sets.newHashSet();
                mJunctionPartitionReadIds.put(unprocessedPartition, readIds);
            }

            readIds.add(readGroup.id());
        }
    }

    private void processCandidateRead(
            final List<String> unprocessedPartitions, final ReadGroup readGroup, final ReadRecord read)
    {
        if(unprocessedPartitions.isEmpty())
            return;

        if(ignoreChromosome(read.Chromosome))
            return;

        /*
        // don't store reads which fall in blacklist regions
        if(!mConfig.RetrieveBlacklistMates && expectedRead != null && mConfig.Blacklist.inBlacklistLocation(
                expectedRead.Chromosome, expectedRead.Position, expectedRead.Position + mConfig.ReadLength))
        {
            return;
        }
        */

        CachedReadGroup cachedReadGroup = null;

        // search all remote partitions for an existing group to add these reads to
        List<String> matchedPartitions = Lists.newArrayList();

        for(String remotePartition : readGroup.remotePartitions())
        {
            Map<String,CachedReadGroup> cachedReadGroups = mCandidatePartitionGroups.get(remotePartition);

            if(cachedReadGroups == null)
                continue;

            if(cachedReadGroup == null)
            {
                cachedReadGroup = cachedReadGroups.get(readGroup.id());

                if(cachedReadGroup != null)
                    matchedPartitions.add(remotePartition);
            }
            else
            {
                if(cachedReadGroups.containsKey(readGroup.id()))
                    matchedPartitions.add(remotePartition);
            }
        }

        if(cachedReadGroup == null)
        {
            cachedReadGroup = new CachedReadGroup();
        }

        cachedReadGroup.Reads.add(read);
        cachedReadGroup.Partitions.addAll(readGroup.remotePartitions());

        // finally ensure each remote unprocessed partition has a link to this cached group
        for(String remotePartition : readGroup.remotePartitions())
        {
            if(!unprocessedPartitions.contains(remotePartition))
                continue;

            if(matchedPartitions.contains(remotePartition))
                continue;

            Map<String,CachedReadGroup> cachedReadGroups = mCandidatePartitionGroups.get(remotePartition);

            if(cachedReadGroups == null)
            {
                cachedReadGroups = Maps.newHashMap();
                mCandidatePartitionGroups.put(remotePartition, cachedReadGroups);
            }

            cachedReadGroups.put(readGroup.id(), cachedReadGroup);
        }
    }

    private void purgePartition(final String chrPartition)
    {
        mJunctionPartitionReadIds.remove(chrPartition); // no further value

        Map<String,CachedReadGroup> cachedReadGroups = mCandidatePartitionGroups.get(chrPartition);

        if(cachedReadGroups == null)
            return;

        mCandidatePartitionGroups.remove(chrPartition);

        // purge any group without unprocessed partitions
        Set<String> purgedGroupReadIds = Sets.newHashSet();
        for(CachedReadGroup cachedReadGroup : cachedReadGroups.values())
        {
            if(cachedReadGroup.Partitions.stream().noneMatch(x -> !mProcessedPartitions.contains(x)))
                purgedGroupReadIds.add(cachedReadGroup.id());
        }

        purgedGroupReadIds.forEach(x -> cachedReadGroups.remove(x));
        mPurgedCandidates += purgedGroupReadIds.size();
    }

    private void logCacheCount(boolean forceLog)
    {
        // int newCount = getCachedReadsCount(null);

        // read groups spanning multiple partitions will be double-counted, but ignore this
        int newCount = mCandidatePartitionGroups.values().stream()
                .mapToInt(x -> x.values().stream().mapToInt(y -> y.Reads.size()).sum()).sum();

        // int newCount = mChrPartitionReadGroupReads.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();

        if(abs(newCount - mLastSnapshotCount) > LOG_CACH_DIFF || forceLog)
        {
            int junctionReadIds = mJunctionPartitionReadIds.values().stream().mapToInt(x -> x.size()).sum();

            SV_LOGGER.info("spanning cache partition processed({}) candidates cached({} -> {} matched={} purged={}) junctionIds({})",
                    mProcessedPartitions.size(), mLastSnapshotCount, newCount, mMatchedCandidates, mPurgedCandidates, junctionReadIds);

            mLastSnapshotCount = newCount;
        }
    }

    private boolean ignoreChromosome(final String chrPartition)
    {
        if(mConfig.SpecificChromosomes.isEmpty())
            return false;

        String chromosome = chrFromChrPartition(chrPartition);
        return !mConfig.SpecificChromosomes.contains(chromosome);
    }

    public synchronized void logStats()
    {
        if(!mConfig.writeReads())
            return;

        logCacheCount(true);

        /*
        int readGroups = 0;

        for(Map<String,List<ReadRecord>> readGroupReads : mChrPartitionReadGroupReads.values())
        {
            readGroups += readGroupReads.size();
        }

        if(mProcessedPartitions.size() > 10)
        {
            SV_LOGGER.info("final spanning cache: readGroups({}) matched({}) purgedCandidate({})",
                    readGroups, mMatchedCandidates, mPurgedCandidates);

            mPerfCounter.logStats();
        }
        */
    }

    private class CachedReadGroup
    {
        public final List<ReadRecord> Reads;
        public final Set<String> Partitions;

        public CachedReadGroup()
        {
            Reads = Lists.newArrayList();
            Partitions = Sets.newHashSet();
        }

        public String id() { return !Reads.isEmpty() ? Reads.get(0).id() : ""; }
        public String toString()
        {
            return format("reads(%s) partitions(%s) id(%s)", Reads.size(), id(), Partitions);
        }
    }
    @VisibleForTesting
    public Map<String,Map<String,List<ReadRecord>>> chrPartitionReadGroupsMap() { return mChrPartitionReadGroupReads; }
    public Map<String,Map<String,CachedReadGroup>> candidatePartitionGroupsMap() { return mCandidatePartitionGroups; }
    public Map<String,Set<String>> junctionPartitionReadIdsMap() { return mJunctionPartitionReadIds; }

    public void reset()
    {
        mChrPartitionReadGroupReads.clear();
        mCandidatePartitionGroups.clear();
        mJunctionPartitionReadIds.clear();
        mProcessedPartitions.clear();
    }

    public int getCachedReadsCount(final String readId)
    {
        Set<CachedReadGroup> cachedReadGroups = Sets.newHashSet();

        mCandidatePartitionGroups.values().stream()
                .filter(x -> readId == null || x.containsKey(readId))
                .forEach(x -> x.values().forEach(y -> cachedReadGroups.add(y)));

        return cachedReadGroups.stream().mapToInt(x -> x.Reads.size()).sum();
    }

    private void processJunctionReadOld(
            final List<String> unprocessedPartitions, final ReadGroup readGroup)
    {
        final ReadRecord read = readGroup.reads().get(0);

        if(ignoreChromosome(read.Chromosome))
            return;

        for(String remotePartition : readGroup.remotePartitions())
        {
            Map<String,List<ReadRecord>> readGroupReads = mChrPartitionReadGroupReads.get(remotePartition);

            if(readGroupReads == null)
                continue;

            List<ReadRecord> cachedReads = readGroupReads.get(readGroup.id());

            if(cachedReads != null)
            {
                cachedReads.forEach(x -> readGroup.addRead(x));
                readGroupReads.remove(readGroup.id());
            }
        }

        // only store this junction readId if other remote unprocessed partitions within this same read group are expected
        if(unprocessedPartitions.isEmpty())
            return;

        // store the junction group's readId against each unprocessed partition, to use to capture expected reads
        for(String unprocessedPartition : unprocessedPartitions)
        {
            Set<String> readIds = mJunctionPartitionReadIds.get(unprocessedPartition);

            if(readIds == null)
            {
                readIds = Sets.newHashSet();
                mJunctionPartitionReadIds.put(unprocessedPartition, readIds);
            }

            readIds.add(readGroup.id());
        }
    }

    private void processCandidateReadOld(
            final List<String> unprocessedPartitions, final ReadGroup readGroup, final ReadRecord read)
    {
        if(unprocessedPartitions.isEmpty())
            return;

        if(ignoreChromosome(read.Chromosome))
            return;

        String readChrPartition = chrPartition(read.Chromosome, read.start());

        /*
        // don't store reads which fall in blacklist regions
        if(!mConfig.RetrieveBlacklistMates && expectedRead != null && mConfig.Blacklist.inBlacklistLocation(
                expectedRead.Chromosome, expectedRead.Position, expectedRead.Position + mConfig.ReadLength))
        {
            return;
        }
        */

        Map<String,List<ReadRecord>> partitionReadGroupReads = mChrPartitionReadGroupReads.get(readChrPartition);

        if(partitionReadGroupReads == null)
        {
            partitionReadGroupReads = Maps.newHashMap();
            mChrPartitionReadGroupReads.put(readChrPartition, partitionReadGroupReads);
        }

        List<ReadRecord> readGroupReads = partitionReadGroupReads.get(readGroup.id());

        if(readGroupReads == null)
        {
            readGroupReads = Lists.newArrayList();
            partitionReadGroupReads.put(readGroup.id(), readGroupReads);
        }

        readGroupReads.add(read);

        // also make a link to this readId from each other expected remote partition back to this read's partition
        for(String unprocessedPartition : unprocessedPartitions)
        {
            Map<String,Set<String>> sourcePartitionMap = mExpectedChrPartitionReadIds.get(unprocessedPartition);
            if(sourcePartitionMap == null)
            {
                sourcePartitionMap = Maps.newHashMap();
                mExpectedChrPartitionReadIds.put(unprocessedPartition, sourcePartitionMap);
            }

            Set<String> readIds = sourcePartitionMap.get(readChrPartition);
            if(readIds == null)
            {
                readIds = Sets.newHashSet();
                sourcePartitionMap.put(readChrPartition, readIds);
            }

            readIds.add(readGroup.id());
        }
    }

    private void purgeUnmatchedReadsOld(final String chrPartition)
    {
        mJunctionPartitionReadIds.remove(chrPartition); // no further value

        Map<String,Set<String>> sourcePartitionMap = mExpectedChrPartitionReadIds.get(chrPartition);
        if(sourcePartitionMap == null)
            return;

        mExpectedChrPartitionReadIds.remove(chrPartition);

        for(Map.Entry<String,Set<String>> entry : sourcePartitionMap.entrySet())
        {
            String sourceChrPartition = entry.getKey();
            Set<String> readIds = entry.getValue();

            Map<String,List<ReadRecord>> readGroupReads = mChrPartitionReadGroupReads.get(sourceChrPartition);

            for(String readId : readIds)
            {
                List<ReadRecord> reads = readGroupReads.get(readId);
                if(reads == null)
                    continue;

                mPurgedCandidates += reads.size();
                readGroupReads.remove(readId);
            }
        }
    }

}
