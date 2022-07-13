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

    private final PerformanceCounter mPerfCounter;
    private final Set<String> mProcessedPartitions;
    private int mLastSnapshotCount;

    public CombinedReadGroups(final SvConfig config)
    {
        mConfig = config;
        mPartitionSize = config.PartitionSize;
        mChrPartitionReadGroupReads = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
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

    private String chrPartition(final String chromosome, int position) { return formChromosomePartition(chromosome, position, mPartitionSize); }

    public synchronized void processSpanningReadGroups(final ChrBaseRegion partitionRegion, final Map<String,ReadGroup> spanningGroups)
    {
        String sourceChrPartition = chrPartition(partitionRegion.Chromosome, partitionRegion.start());
        mProcessedPartitions.add(sourceChrPartition);

        Set<ExpectedRead> addedReads = Sets.newHashSet();

        // look for reads which have already been found (and therefore written)
        for(Map.Entry<String,ReadGroup> entry : spanningGroups.entrySet())
        {
            ReadGroup readGroup = entry.getValue();

            for(ReadRecord read : readGroup.reads())
            {
                processNewRead(sourceChrPartition, readGroup, read, null, addedReads);
            }

            // also the expected reads
            List<ExpectedRead> expectedReads = ExpectedRead.getExpectedReads(readGroup);

            for(ExpectedRead read : expectedReads)
            {
                processNewRead(sourceChrPartition, readGroup, null, read, addedReads);
            }
        }

        // purge any expected reads which were already found
        purgeUnmatchedReads(sourceChrPartition, addedReads);

        setCacheCount();
    }

    private void processNewRead(
            final String sourceChrPartition, final ReadGroup readGroup, final ReadRecord actualRead, final ExpectedRead expectedRead,
            final Set<ExpectedRead> addedReads)
    {
        String readChromosome = actualRead != null ? actualRead.Chromosome : expectedRead.Chromosome;
        int readPosition = actualRead != null ? actualRead.start() : expectedRead.Position;

        if(ignoreChromosome(readChromosome))
            return;

        String readChrPartition = chrPartition(readChromosome, readPosition);
        ExpectedRead read = expectedRead != null ? expectedRead : ExpectedRead.fromRead(actualRead);

        if(findExistingRead(readChrPartition, readGroup, actualRead, read))
            return;

        // only store this read against partitions other than the one it came through with, and factor in remote partitions already processed
        int unprocessedPartitions = 0;
        for(String remotePartition : readGroup.remotePartitions())
        {
            String remoteChr = chromosomeFromChromosomePartition(remotePartition);

            if(ignoreChromosome(remoteChr))
                continue;

            if(mProcessedPartitions.contains(remotePartition))
                continue;

            ++unprocessedPartitions;
            Map<String,List<ExpectedRead>> partitionReadGroupReads = mChrPartitionReadGroupReads.get(remotePartition);

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
            read.setExpectedMatchCount(unprocessedPartitions);
            readGroupReads.add(read);
            addedReads.add(read);
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

        return true;
    }

    private void purgeUnmatchedReads(final String chrPartition, final Set<ExpectedRead> addedReads)
    {
        Map<String,List<ExpectedRead>> readGroupReads = mChrPartitionReadGroupReads.get(chrPartition);

        if(readGroupReads == null)
            return;

        Set<String> emptyGroups = Sets.newHashSet();
        for(Map.Entry<String,List<ExpectedRead>> entry : readGroupReads.entrySet())
        {
            List<ExpectedRead> reads = entry.getValue();

            int index = 0;
            while(index < reads.size())
            {
                ExpectedRead read = reads.get(index);

                if(read.Found && !addedReads.contains(read))
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
                emptyGroups.add(entry.getKey());
        }

        emptyGroups.forEach(x -> readGroupReads.remove(x));

        if(readGroupReads.isEmpty())
        {
            mChrPartitionReadGroupReads.remove(chrPartition);
            return;
        }
    }

    private void setCacheCount()
    {
        int newCount = mChrPartitionReadGroupReads.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();

        if(abs(newCount - mLastSnapshotCount) > LOG_CACH_DIFF)
        {
            mLastSnapshotCount = newCount;

            SV_LOGGER.info("spanning partition processed({}) groups cached({} -> {})",
                    mProcessedPartitions.size(), mLastSnapshotCount, newCount);
        }
    }

    private boolean ignoreChromosome(final String chrPartition)
    {
        if(mConfig.SpecificChromosomes.isEmpty())
            return false;

        String chromosome = chromosomeFromChromosomePartition(chrPartition);
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
                List<ExpectedRead> foundReads = entry.getValue().stream().filter(x -> x.Found).collect(Collectors.toList());
                List<ExpectedRead> unfoundReads = entry.getValue().stream().filter(x -> !x.Found).collect(Collectors.toList());

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

        SV_LOGGER.info("final spanning partition cache: partitions({}) readGroups({}) reads({} found={} nonSuppMissed={} suppMissed={})",
                partitions, readGroups, totalCachedReads, foundReadCount, missedNonSuppReadCount, missedSuppReadCount);

        mPerfCounter.logStats();
    }
}
