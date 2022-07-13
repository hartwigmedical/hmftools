package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.BAM;
import static com.hartwig.hmftools.svprep.WriteType.READS;

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

public class CombinedReadGroups
{
    private final SvConfig mConfig;
    private final int mPartitionSize;

    private final Map<String,Map<String,List<ExpectedRead>>> mChrPartitionReadGroupReads; // keyed by chromosome-partition then readId

    // a map of remote chr-partition to read partition to readIds
    private final Map<String,Map<String,Set<String>>> mExpectedChrPartitionReadIds;

    private final PerformanceCounter mPerfCounter;
    private final Set<String> mProcessedPartitions;

    private int mLastSnapshotCount;
    private int mMatchedGroups;

    // logic to keep expected reads for a final search for them in the BAM has been disabled

    public CombinedReadGroups(final SvConfig config)
    {
        mConfig = config;
        mPartitionSize = config.PartitionSize;
        mChrPartitionReadGroupReads = Maps.newHashMap();
        mExpectedChrPartitionReadIds = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mLastSnapshotCount = 0;
        mMatchedGroups = 0;
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

    public synchronized void processSpanningReadGroups(
            final ChrBaseRegion partitionRegion, final Map<String,ReadGroup> spanningGroups,
            final Map<String,List<ExpectedRead>> missedReadsMap)
    {
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

        // only store this read if other remote partitions are expected and factor in remote partitions already processed
        if(unprocessedPartitions.isEmpty())
            return;

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

        // also make a link to this readId from the expected remote partition back to this read's partition
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

                    if(!read.found())
                    {
                        // collect missed reads to then query the BAM for them
                        List<ExpectedRead> missedReads = missedReadsMap.get(readId);
                        if(missedReads == null)
                        {
                            missedReads = Lists.newArrayList();
                            missedReadsMap.put(readId, missedReads);
                        }
                        missedReads.add(read);
                    }

                    if(!addedReads.contains(read)) // read.found() &&
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
            SV_LOGGER.info("spanning partition processed({}) groups cached({} -> {}) matchedGroups({})",
                    mProcessedPartitions.size(), mLastSnapshotCount, newCount, mMatchedGroups);

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

    public void writeRemainingReadGroups(final ResultsWriter writer, final Set<WriteType> writeTypes)
    {
        if(!writeTypes.contains(BAM) && !writeTypes.contains(READS))
            return;

        int readGroups = 0;
        int foundReadCount = 0;
        int missedSuppReadCount = 0;
        int missedNonSuppReadCount = 0;

        Map<String,List<ExpectedRead>> foundReadGroups = Maps.newHashMap();
        Map<String,List<ExpectedRead>> missedReadGroups = Maps.newHashMap();

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

                if(!foundReads.isEmpty())
                    foundReadGroups.put(entry.getKey(), foundReads);

                if(!unfoundReads.isEmpty())
                    missedReadGroups.put(entry.getKey(), unfoundReads);
            }
        }

        int partitions = mChrPartitionReadGroupReads.size();
        int totalCachedReads = foundReadCount + missedNonSuppReadCount + missedSuppReadCount;
        int totalMissed = missedNonSuppReadCount + missedSuppReadCount;

        SV_LOGGER.info("final spanning partition cache: partitions({}) readGroups({}) matched({}) reads({} found={} missed={})",
                partitions, readGroups, mMatchedGroups, totalCachedReads, foundReadCount, totalMissed);

        mPerfCounter.logStats();
    }

    @VisibleForTesting
    public Map<String,Map<String,List<ExpectedRead>>> chrPartitionReadGroupsMap() { return mChrPartitionReadGroupReads; }

    public void reset()
    {
        mChrPartitionReadGroupReads.clear();
        mProcessedPartitions.clear();
    }
}
