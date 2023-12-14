package com.hartwig.hmftools.markdups;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.APP_NAME;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.addConfig;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.markdups.common.Constants.LOCK_ACQUIRE_LONG_TIME_MS;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSampler;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.common.Statistics;
import com.hartwig.hmftools.markdups.consensus.ConsensusReads;
import com.hartwig.hmftools.markdups.write.BamWriter;
import com.hartwig.hmftools.markdups.write.FileWriterCache;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

public class MarkDuplicates
{
    private final MarkDupsConfig mConfig;

    public MarkDuplicates(final ConfigBuilder configBuilder)
    {
        mConfig = new MarkDupsConfig(configBuilder);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        MD_LOGGER.info("sample({}) starting mark duplicates", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        setReadLength();

        FileWriterCache fileWriterCache = new FileWriterCache(mConfig);

        PartitionDataStore partitionDataStore = new PartitionDataStore(mConfig);

        // partition all chromosomes
        Queue<ChrBaseRegion> partitions = new ConcurrentLinkedQueue<>();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            List<ChrBaseRegion> chrPartitions = partitionChromosome(
                    chromosomeStr, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize);

            partitions.addAll(chrPartitions);
        }

        List<PartitionThread> partitionTasks = Lists.newArrayList();
        List<Thread> workers = new ArrayList<>();

        int partitionCount = partitions.size();

        for(int i = 0; i < min(partitionCount, mConfig.Threads); ++i)
        {
            PartitionThread partitionThread = new PartitionThread(i, mConfig, partitions, fileWriterCache, partitionDataStore);
            partitionTasks.add(partitionThread);
            workers.add(partitionThread);
        }

        MD_LOGGER.debug("splitting {} partitions across {} threads", partitionCount, partitionTasks.size());

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                MD_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        MD_LOGGER.info("all partition tasks complete");

        long unmappedReads = writeUnmappedReads(fileWriterCache);

        List<PartitionReader> partitionReaders = partitionTasks.stream().map(x -> x.partitionReader()).collect(Collectors.toList());

        int maxLogFragments = (mConfig.RunChecks || mConfig.LogFinalCache) ? 100 : 0;
        int totalUnwrittenFragments = 0;
        ConsensusReads consensusReads = new ConsensusReads(mConfig.RefGenome);
        consensusReads.setDebugOptions(mConfig.RunChecks);

        // write any orphaned or remaining fragments (can be supplementaries)
        BamWriter recordWriter = partitionTasks.get(0).bamWriter();

        for(PartitionData partitionData : partitionDataStore.partitions())
        {
            int cachedReadCount = partitionData.writeRemainingReads(recordWriter, consensusReads, maxLogFragments > 0);
            totalUnwrittenFragments += cachedReadCount;
            maxLogFragments = max(0, maxLogFragments - cachedReadCount);
        }

        if(totalUnwrittenFragments > 0)
        {
            MD_LOGGER.info("wrote {} remaining cached fragments", totalUnwrittenFragments);
        }

        List<PerformanceCounter> combinedPerfCounters = mergePerfCounters(partitionReaders);

        Statistics combinedStats = new Statistics();
        partitionReaders.forEach(x -> combinedStats.merge(x.statistics()));
        partitionDataStore.partitions().forEach(x -> combinedStats.merge(x.statistics()));
        combinedStats.ConsensusStats.merge(consensusReads.consensusStats());

        // free up any processing state
        partitionReaders.clear();

        fileWriterCache.close();

        if(fileWriterCache.runSortMergeIndex())
        {
            // log interim time
            MD_LOGGER.info("BAM duplicate processing complete, mins({})", runTimeMinsStr(startTimeMs));

            if(!fileWriterCache.sortAndIndexBams())
            {
                MD_LOGGER.error("sort-merge-index failed");
                System.exit(1);
            }
        }

        combinedStats.logStats();

        long totalWrittenReads = fileWriterCache.totalWrittenReads();
        long unmappedDroppedReads = mConfig.UnmapRegions.stats().SupplementaryCount.get();

        if(mConfig.UnmapRegions.enabled())
        {
            MD_LOGGER.info("unmapped stats: {}", mConfig.UnmapRegions.stats().toString());
        }

        if(combinedStats.TotalReads + unmappedReads != totalWrittenReads + unmappedDroppedReads)
        {
            long difference = combinedStats.TotalReads + unmappedReads - totalWrittenReads - unmappedDroppedReads;

            MD_LOGGER.warn("reads processed({}) vs written({}) mismatch diffLessDropped({})",
                    combinedStats.TotalReads + unmappedReads, totalWrittenReads, difference);
        }

        if(mConfig.WriteStats)
        {
            combinedStats.writeDuplicateStats(mConfig);

            if(mConfig.UMIs.Enabled)
            {
                combinedStats.UmiStats.writePositionFragmentsData(mConfig);

                if(mConfig.UMIs.BaseStats)
                {
                    combinedStats.UmiStats.writeUmiBaseDiffStats(mConfig);
                    combinedStats.UmiStats.writeUmiBaseFrequencyStats(mConfig);
                }
            }
        }

        logPerformanceStats(combinedPerfCounters, partitionDataStore);

        MD_LOGGER.info("Mark duplicates complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private long writeUnmappedReads(final FileWriterCache fileWriterCache)
    {
        if(mConfig.SpecificChrRegions.hasFilters() || !mConfig.WriteBam)
            return 0;

        BamWriter bamWriter = fileWriterCache.getFullyUnmappedReadsBamWriter();

        BamReader bamReader = new BamReader(mConfig);

        AtomicLong unmappedCount = new AtomicLong();

        bamReader.queryUnmappedReads((final SAMRecord record) ->
        {
            bamWriter.writeRead(record, FragmentStatus.UNSET);
            unmappedCount.incrementAndGet();
        });

        if(unmappedCount.get() > 0)
        {
            MD_LOGGER.debug("wrote {} unmapped reads", unmappedCount);
        }

        return unmappedCount.get();
    }

    private void setReadLength()
    {
        if(mConfig.readLength() > 0) // skip if set in config
            return;

        // sample the BAM to determine read length
        BamSampler bamSampler = new BamSampler(mConfig.RefGenomeFile);

        ChrBaseRegion sampleRegion = !mConfig.SpecificChrRegions.Regions.isEmpty() ?
                mConfig.SpecificChrRegions.Regions.get(0) : bamSampler.defaultRegion();

        int readLength = DEFAULT_READ_LENGTH;

        if(bamSampler.calcBamCharacteristics(mConfig.BamFiles.get(0), sampleRegion) && bamSampler.maxReadLength() > 0)
        {
            readLength = bamSampler.maxReadLength();
            MD_LOGGER.info("BAM sampled max read-length({})", readLength);
        }
        else
        {
            MD_LOGGER.warn("BAM read-length sampling failed, using default read length({})", DEFAULT_READ_LENGTH);
        }

        mConfig.setReadLength(readLength);
    }

    private List<PerformanceCounter> mergePerfCounters(final List<PartitionReader> partitionReaders)
    {
        List<PerformanceCounter> combinedPerfCounters = partitionReaders.get(0).perfCounters();

        for(int i = 1; i < partitionReaders.size(); ++i)
        {
            List<PerformanceCounter> chrPerfCounters = partitionReaders.get(i).perfCounters();

            for(int j = 0; j < chrPerfCounters.size(); ++j)
            {
                combinedPerfCounters.get(j).merge(chrPerfCounters.get(j));
            }
        }

        return combinedPerfCounters;
    }

    private void logPerformanceStats(final List<PerformanceCounter> combinedPerfCounters, final PartitionDataStore partitionDataStore)
    {
        if(mConfig.PerfDebug)
        {
            for(int j = 0; j < combinedPerfCounters.size(); ++j)
            {
                PerformanceCounter perfCounter = combinedPerfCounters.get(j);

                if(j == 0)
                    perfCounter.logIntervalStats(10);
                else
                    perfCounter.logStats();
            }

            // check partition store locking times
            double totalLockTimeMs = 0;

            for(PartitionData partitionData : partitionDataStore.partitions())
            {
                double lockTime = partitionData.totalLockTimeMs();

                totalLockTimeMs += lockTime;

                if(lockTime > LOCK_ACQUIRE_LONG_TIME_MS)
                {
                    MD_LOGGER.debug("partition({}) lock-acquisition time({}ms)",
                            partitionData.partitionStr(), format("%.1f", lockTime));
                }
            }

            if(totalLockTimeMs > LOCK_ACQUIRE_LONG_TIME_MS)
            {
                MD_LOGGER.debug("partition cache total lock-acquisition time({}s)",
                        format("%.3f", totalLockTimeMs / 1000));
            }
        }
        else
        {
            combinedPerfCounters.forEach(x -> x.logStats());
        }

    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        MarkDuplicates markDuplicates = new MarkDuplicates(configBuilder);
        markDuplicates.run();
    }
}
