package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.EXPECTED;

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
import com.hartwig.hmftools.svprep.reads.ExpectedRead;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.ReadType;

public class CombinedReadGroups
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

    private final Map<String,Map<String,List<ExpectedRead>>> mChrPartitionReadGroupReads; // keyed by chromosome-partition then readId

    // a map of remote chr-partition to read partition to readIds
    private final Map<String,Map<String,Set<String>>> mExpectedChrPartitionReadIds;

    private final Set<String> mUnmappedReadIds;

    private final PerformanceCounter mPerfCounter;
    private final Set<String> mProcessedPartitions;

    private int mLastSnapshotCount;
    private int mMatchedGroups;
    private int mPurgedCandidates;

    public CombinedReadGroups(final SvConfig config)
    {
        mConfig = config;
        mPartitionSize = config.PartitionSize;
        mChrPartitionReadGroupReads = Maps.newHashMap();
        mExpectedChrPartitionReadIds = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mUnmappedReadIds = Sets.newHashSet();
        mLastSnapshotCount = 0;
        mMatchedGroups = 0;
        mPurgedCandidates = 0;
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

    private static String chrFromChrPartition(final String chrPartition)
    {
        return chrPartition.split(CHR_PARTITION_DELIM, 2)[0];
    }

    private String chrPartition(final String chromosome, int position) { return formChromosomePartition(chromosome, position, mPartitionSize); }

    public synchronized Set<String> getExpectedReadIds(final ChrBaseRegion partitionRegion)
    {
        Set<String> expectedReadIds = Sets.newHashSet();

        String chrPartition = chrPartition(partitionRegion.Chromosome, partitionRegion.start());
        Map<String,Set<String>> sourcePartitionMap = mExpectedChrPartitionReadIds.get(chrPartition);
        if(sourcePartitionMap == null)
            return expectedReadIds;

        // only include readIds from groups with a written read, ie not remote candidates
        for(Map.Entry<String,Set<String>> entry : sourcePartitionMap.entrySet())
        {
            String sourceChrPartition = entry.getKey();
            Set<String> readIds = entry.getValue();

            Map<String,List<ExpectedRead>> readGroupReads = mChrPartitionReadGroupReads.get(sourceChrPartition);

            for(String readId : readIds)
            {
                List<ExpectedRead> reads = readGroupReads.get(readId);
                if(reads == null || reads.stream().allMatch(x -> x.remoteCandidateGroup()))
                    continue;

                expectedReadIds.add(readId);
            }
        }

        // sourcePartitionMap.values().forEach(x -> expectedReadIds.addAll(x));

        return expectedReadIds;
    }

    public synchronized void processSpanningReadGroups(
            final ChrBaseRegion partitionRegion, final Map<String,ReadGroup> spanningGroups,
            final Map<String,List<ExpectedRead>> missedReadsMap)
    {
        mPerfCounter.start();

        String sourceChrPartition = chrPartition(partitionRegion.Chromosome, partitionRegion.start());
        mProcessedPartitions.add(sourceChrPartition);

        Set<ExpectedRead> addedReads = Sets.newHashSet();

        // look for reads which have already been found (and therefore written)
        for(Map.Entry<String,ReadGroup> entry : spanningGroups.entrySet())
        {
            ReadGroup readGroup = entry.getValue();

            List<String> unprocessedPartitions = Lists.newArrayList();
            for(String remotePartition : readGroup.remotePartitions())
            {
                if(mProcessedPartitions.contains(remotePartition))
                    continue;

                if(!mConfig.SpecificChromosomes.isEmpty() && ignoreChromosome(chrFromChrPartition(remotePartition)))
                    continue;

                unprocessedPartitions.add(remotePartition);
            }

            for(ReadRecord read : readGroup.reads())
            {
                processNewRead(unprocessedPartitions, readGroup, read, null, addedReads);
            }

            // also the expected reads
            List<ExpectedRead> expectedReads = ExpectedRead.getExpectedReads(readGroup);

            for(ExpectedRead read : expectedReads)
            {
                processNewRead(unprocessedPartitions, readGroup, null, read, addedReads);
            }
        }

        // purge any expected reads which were already found
        purgeUnmatchedReads(sourceChrPartition, addedReads, missedReadsMap);

        setCacheCount();

        mPerfCounter.stop();
    }

    private void processNewRead(
            final List<String> unprocessedPartitions, final ReadGroup readGroup, final ReadRecord actualRead,
            final ExpectedRead expectedRead, final Set<ExpectedRead> addedReads)
    {
        String readChromosome = actualRead != null ? actualRead.Chromosome : expectedRead.Chromosome;
        int readPosition = actualRead != null ? actualRead.start() : expectedRead.Position;

        if(ignoreChromosome(readChromosome))
            return;

        String readChrPartition = chrPartition(readChromosome, readPosition);
        ExpectedRead read = expectedRead != null ? expectedRead : ExpectedRead.fromRead(actualRead);

        if(findExistingRead(readChrPartition, readGroup, actualRead, read))
            return;

        // supplementaries logic:
        // an expected supplementary read, if matched, will collect its actual read that was cached to then be written to the BAM
        // an actual supplementary read, if matched, is marked as belonging to a read group with non-supplementaries
        // an actual supplementary read, if not matched, caches its actual read with the stored expected read
        // unclaimed supplementaries are not given back in the purge process since their non-supps are likely duplicates

        // only store this read if other remote partitions within this same read group are expected
        // and factor in remote partitions already processed
        if(unprocessedPartitions.isEmpty())
            return;

        // don't store reads which fall in blacklist regions
        if(!mConfig.RetrieveBlacklistMates && expectedRead != null && mConfig.Blacklist.inBlacklistLocation(
                expectedRead.Chromosome, expectedRead.Position, expectedRead.Position + mConfig.ReadLength))
        {
            return;
        }

        Map<String,List<ExpectedRead>> partitionReadGroupReads = mChrPartitionReadGroupReads.get(readChrPartition);

        if(partitionReadGroupReads == null)
        {
            partitionReadGroupReads = Maps.newHashMap();
            mChrPartitionReadGroupReads.put(readChrPartition, partitionReadGroupReads);
        }

        List<ExpectedRead> readGroupReads = partitionReadGroupReads.get(readGroup.id());

        if(readGroupReads == null)
        {
            readGroupReads = Lists.newArrayList();
            partitionReadGroupReads.put(readGroup.id(), readGroupReads);
        }

        // add this new read
        read.setExpectedMatchCount(unprocessedPartitions.size());
        readGroupReads.add(read);
        addedReads.add(read);

        // for groups which are conditional to be written on their mate being a junction or support read, cache their read
        // this also prevents supplementaries of duplicates being written
        if(readGroup.conditionalOnRemoteReads())
        {
            if(actualRead != null)
                read.setCachedRead(actualRead);

            read.markRemoteCandidateGroup();
            // return; // don't register the read ID as expected
        }

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

    private boolean findExistingRead(
            final String readChrPartition, final ReadGroup readGroup, final ReadRecord actualRead, final ExpectedRead read)
    {
        Map<String,List<ExpectedRead>> readGroupReads = mChrPartitionReadGroupReads.get(readChrPartition);

        if(readGroupReads == null)
            return false;

        List<ExpectedRead> existingGroupReads = readGroupReads.get(readGroup.id());

        if(existingGroupReads == null)
            return false;

        ExpectedRead matchedRead = existingGroupReads.stream().filter(x -> x.matches(read)).findFirst().orElse(null);
        if(matchedRead == null)
            return false;

        if(actualRead != null)
        {
            if(matchedRead.found())
                actualRead.setWritten();
            else
                matchedRead.markFound();

            if(readGroup.conditionalOnRemoteReads() && !matchedRead.remoteCandidateGroup())
                readGroup.markHasRemoteJunctionReads();
            else if(matchedRead.remoteCandidateGroup() && actualRead.readType() == ReadType.EXPECTED)
                readGroup.markConditionalOnRemoteReads();
        }
        else
        {
            if(matchedRead.remoteCandidateGroup() && readGroup.reads().stream().allMatch(x -> x.readType() == EXPECTED))
            {
                readGroup.markConditionalOnRemoteReads();
            }
            else if(matchedRead.hasCachedRead() && !readGroup.conditionalOnRemoteReads())
            {
                // pick up supplementaries or candidate supporting reads which were waiting for their group
                readGroup.addRead(matchedRead.getCachedRead());
                readGroup.setGroupState();
            }
        }

        matchedRead.setExpectedMatchCount(readGroup.partitionCount() - 1);
        matchedRead.registerMatch();

        if(matchedRead.fullyMatched())
        {
            existingGroupReads.remove(matchedRead);

            if(existingGroupReads.isEmpty())
            {
                readGroupReads.remove(readGroup.id());
                ++mMatchedGroups;
            }
        }

        return true;
    }

    private void purgeUnmatchedReads(
            final String chrPartition, final Set<ExpectedRead> addedReads, final Map<String,List<ExpectedRead>> missedReadsMap)
    {
        Map<String,Set<String>> sourcePartitionMap = mExpectedChrPartitionReadIds.get(chrPartition);
        if(sourcePartitionMap == null)
            return;

        mExpectedChrPartitionReadIds.remove(chrPartition); // no further value

        for(Map.Entry<String,Set<String>> entry : sourcePartitionMap.entrySet())
        {
            String sourceChrPartition = entry.getKey();
            Set<String> readIds = entry.getValue();

            Map<String,List<ExpectedRead>> readGroupReads = mChrPartitionReadGroupReads.get(sourceChrPartition);

            Set<String> emptyGroups = Sets.newHashSet();

            for(String readId : readIds)
            {
                List<ExpectedRead> reads = readGroupReads.get(readId);
                if(reads == null)
                    continue;

                int index = 0;
                while(index < reads.size())
                {
                    ExpectedRead read = reads.get(index);

                    if(read.hasCachedRead() || read.remoteCandidateGroup())
                    {
                        // unmatched supplementaries likely belong to duplicates are so are dropped now with retrieval of their non-supps
                        // otherwise are lone candidate reads cached to avoid a re-slice but can now be purged
                        reads.remove(index);
                        ++mPurgedCandidates;
                        continue;
                    }

                    if(!read.found())
                    {
                        // collect missed reads to then query the BAM for them back in the partition slicer
                        List<ExpectedRead> missedReads = missedReadsMap.get(readId);
                        if(missedReads == null)
                        {
                            missedReads = Lists.newArrayList();
                            missedReadsMap.put(readId, missedReads);
                        }
                        missedReads.add(read);
                    }

                    if(!addedReads.contains(read))
                    {
                        read.registerMatch();

                        if(read.fullyMatched())
                        {
                            reads.remove(index);
                            continue;
                        }
                    }

                    ++index;
                }

                if(reads.isEmpty())
                    emptyGroups.add(readId);
            }

            emptyGroups.forEach(x -> readGroupReads.remove(x));
        }
    }

    private void setCacheCount()
    {
        int newCount = mChrPartitionReadGroupReads.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();

        if(abs(newCount - mLastSnapshotCount) > LOG_CACH_DIFF)
        {
            SV_LOGGER.info("spanning partition processed({}) groups cached({} -> {}) matchedGroups({}) purgedCandidates({}) unmapped({})",
                    mProcessedPartitions.size(), mLastSnapshotCount, newCount, mMatchedGroups, mPurgedCandidates, mUnmappedReadIds.size());

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

    public void logStats(final ResultsWriter writer)
    {
        if(!mConfig.writeReads())
            return;

        int readGroups = 0;
        int foundReadCount = 0;
        int missedSuppReadCount = 0;
        int missedNonSuppReadCount = 0;

        for(Map<String,List<ExpectedRead>> readGroupReads : mChrPartitionReadGroupReads.values())
        {
            readGroups += readGroupReads.size();

            for(Map.Entry<String,List<ExpectedRead>> entry : readGroupReads.entrySet())
            {
                List<ExpectedRead> foundReads = entry.getValue().stream().filter(x -> x.found()).collect(Collectors.toList());
                List<ExpectedRead> unfoundReads = entry.getValue().stream().filter(x -> !x.found()).collect(Collectors.toList());

                foundReadCount += foundReads.size();
                missedSuppReadCount += unfoundReads.stream().filter(x -> x.IsSupplementary).count();
                missedNonSuppReadCount += unfoundReads.stream().filter(x -> !x.IsSupplementary).count();
            }
        }

        int totalCachedReads = foundReadCount + missedNonSuppReadCount + missedSuppReadCount;
        int totalMissed = missedNonSuppReadCount + missedSuppReadCount;

        if(mProcessedPartitions.size() > 10 || totalMissed > 0)
        {
            SV_LOGGER.info("final spanning partition cache: readGroups({}) matched({}) purgedCandidate({}) reads({} found={} missed={})",
                    readGroups, mMatchedGroups, mPurgedCandidates, totalCachedReads, foundReadCount, totalMissed);

            mPerfCounter.logStats();
        }
    }

    @VisibleForTesting
    public Map<String,Map<String,List<ExpectedRead>>> chrPartitionReadGroupsMap() { return mChrPartitionReadGroupReads; }

    public void reset()
    {
        mChrPartitionReadGroupReads.clear();
        mProcessedPartitions.clear();
    }
}
