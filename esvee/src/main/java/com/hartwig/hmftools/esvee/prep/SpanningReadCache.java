package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;

public class SpanningReadCache
{
    private final PrepConfig mConfig;
    private final int mPartitionSize;

    // a cache of read groups found from each chr-partition to help retrieve complete fragments (ie expected reads within a group)
    // when a new chr-partition completes, the following steps are done:
    // - record that the chr-partition has been processed (in mProcessedPartitions)
    // for read groups that support a junction:
    // - pick up an cached candidate or supplementary reads, which then will be written to file
    // - for the group's unprocessed remote partitions, cache the readId to aid with identifying expected reads
    // for read groups with only candidate / supplementary reads
    // - if no unprocessed partitions then drop immediately
    // - otherwise cache as now (ie source partition and readId)

    private final Map<String,Map<String, CachedReadGroup>> mCandidatePartitionGroups; // keyed by chromosome-partition then readId
    private final Map<String,Set<String>> mJunctionPartitionReadIds;

    private final PerformanceCounter mPerfCounter;
    private final Set<String> mProcessedPartitions;

    private int mLastSnapshotCount;
    private int mMatchedCandidates;
    private int mPurgedCandidates;
    private final CandidateBamWriter mCandidateBamWriter;

    public SpanningReadCache(final PrepConfig config)
    {
        mConfig = config;
        mPartitionSize = config.PartitionSize;
        mCandidatePartitionGroups = Maps.newHashMap();
        mJunctionPartitionReadIds = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mLastSnapshotCount = 0;
        mMatchedCandidates = 0;
        mPurgedCandidates = 0;
        mCandidateBamWriter = new CandidateBamWriter(config);
        mPerfCounter = new PerformanceCounter("SpanningReads");
    }

    private static final String CHR_PARTITION_DELIM = "_";
    private static final int LOG_CACH_DIFF = 50000;

    public CandidateBamWriter candidateBamWriter() { return mCandidateBamWriter; }

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }

    public static String chrFromChrPartition(final String chrPartition)
    {
        return chrPartition.split(CHR_PARTITION_DELIM, 2)[0];
    }

    private String chrPartition(final String chromosome, int position) { return formChromosomePartition(chromosome, position, mPartitionSize); }

    public synchronized Set<String> getExpectedReadIds(final ChrBaseRegion partitionRegion)
    {
        String chrPartition = chrPartition(partitionRegion.Chromosome, partitionRegion.start());
        Set<String> expectedReadIds = mJunctionPartitionReadIds.get(chrPartition);
        return expectedReadIds != null ? Sets.newHashSet(expectedReadIds) : Sets.newHashSet();
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
                    .filter(x -> mConfig.SpecificChrRegions.Chromosomes.isEmpty() || !ignoreChromosome(chrFromChrPartition(x)))
                    .collect(Collectors.toList());

            if(readGroup.conditionalOnRemoteReads())
            {
                for(PrepRead read : readGroup.reads())
                {
                    processCandidateRead(unprocessedPartitions, readGroup, read);
                }
            }
            else
            {
                processJunctionRead(unprocessedPartitions, readGroup, sourceChrPartition);
            }
        }

        // purge any cached candidate reads and junction readIds which are no longer relevant
        purgePartition(sourceChrPartition);

        logCacheCount(false);

        mPerfCounter.stop();
    }

    private void processJunctionRead(
            final List<String> unprocessedPartitions, final ReadGroup readGroup, final String sourceChrPartition)
    {
        final PrepRead read = readGroup.reads().get(0);

        if(ignoreChromosome(read.Chromosome))
            return;

        if(mConfig.UseCacheBam)
        {
            mCandidateBamWriter.addJunctionReadId(readGroup.remotePartitions(), readGroup.id());
        }
        else
        {
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
            final List<String> unprocessedPartitions, final ReadGroup readGroup, final PrepRead read)
    {
        if(unprocessedPartitions.isEmpty())
            return;

        if(ignoreChromosome(read.Chromosome))
            return;

        if(mConfig.UseCacheBam)
        {
            mCandidateBamWriter.writeCandidateRead(read);
            return;
        }

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
            cachedReadGroup = new CachedReadGroup(readGroup.id());
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
                purgedGroupReadIds.add(cachedReadGroup.ReadId);
        }

        purgedGroupReadIds.forEach(x -> cachedReadGroups.remove(x));
        mPurgedCandidates += purgedGroupReadIds.size();
    }

    private void logCacheCount(boolean forceLog)
    {
        if(mConfig.UseCacheBam)
        {
            if(!forceLog)
                return;

            int junctionReadIds = mJunctionPartitionReadIds.values().stream().mapToInt(x -> x.size()).sum();

            SV_LOGGER.info("spanning cache partition processed({}) junctionIds({})", mProcessedPartitions.size(), junctionReadIds);
            return;
        }

        // read groups spanning multiple partitions will be double-counted, but ignore this
        int newCount = mCandidatePartitionGroups.values().stream()
                .mapToInt(x -> x.values().stream().mapToInt(y -> y.Reads.size()).sum()).sum();

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
        if(mConfig.SpecificChrRegions.Chromosomes.isEmpty())
            return false;

        String chromosome = chrFromChrPartition(chrPartition);
        return !mConfig.SpecificChrRegions.Chromosomes.contains(chromosome);
    }

    public synchronized void logStats()
    {
        if(!mConfig.PerfDebug)
            return;

        logCacheCount(true);
        mPerfCounter.logStats();
    }

    private class CachedReadGroup
    {
        public final String ReadId;
        public final List<PrepRead> Reads;
        public final Set<String> Partitions;

        public CachedReadGroup(final String readId)
        {
            ReadId = readId;
            Reads = Lists.newArrayList();
            Partitions = Sets.newHashSet();
        }

        public String toString()
        {
            return format("reads(%s) partitions(%s) id(%s)", Reads.size(), ReadId, Partitions);
        }
    }

    @VisibleForTesting
    public Map<String,Set<String>> junctionPartitionReadIdsMap() { return mJunctionPartitionReadIds; }

    public void reset()
    {
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
}
