package com.hartwig.hmftools.markdups;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.APP_NAME;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.addConfig;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.common.Statistics;
import com.hartwig.hmftools.markdups.consensus.ConsensusReads;

import org.jetbrains.annotations.NotNull;

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

        List<ChromosomeReader> chromosomeReaders = Lists.newArrayList();

        RefGenomeCoordinates refGenomeCoordinates = mConfig.RefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        FileWriterCache fileWriterCache = new FileWriterCache(mConfig);

        PartitionDataStore partitionDataStore = new PartitionDataStore(mConfig);
        final List<Callable> callableList = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            ChrBaseRegion chrBaseRegion = new ChrBaseRegion(chromosomeStr, 1, refGenomeCoordinates.Lengths.get(chromosome));

            ChromosomeReader chromosomeReader = new ChromosomeReader(chrBaseRegion, mConfig, fileWriterCache, partitionDataStore);
            chromosomeReaders.add(chromosomeReader);
            callableList.add(chromosomeReader);
        }

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        int maxLogFragments = (mConfig.RunChecks || mConfig.LogFinalCache) ? 100 : 0;
        int totalUnwrittenFragments = 0;
        ConsensusReads consensusReads = new ConsensusReads(mConfig.RefGenome);
        consensusReads.setDebugOptions(mConfig.RunChecks);

        BamWriter recordWriter = chromosomeReaders.get(0).recordWriter();

        for(PartitionData partitionData : partitionDataStore.partitions())
        {
            int cachedReadCount = partitionData.writeRemainingReads(recordWriter, consensusReads, maxLogFragments > 0);
            totalUnwrittenFragments += cachedReadCount;
            maxLogFragments = max(0, maxLogFragments - cachedReadCount);
        }

        if(totalUnwrittenFragments > 0)
        {
            MD_LOGGER.info("wrote {} cached fragments", totalUnwrittenFragments);
        }

        fileWriterCache.close();

        MD_LOGGER.debug("all chromosome tasks complete");

        Statistics combinedStats = new Statistics();
        chromosomeReaders.forEach(x -> combinedStats.merge(x.statistics()));
        partitionDataStore.partitions().forEach(x -> combinedStats.merge(x.statistics()));

        int totalWrittenReads = chromosomeReaders.stream().mapToInt(x -> x.recordWriter().recordWriteCount()).sum();

        combinedStats.logStats();

        if(combinedStats.TotalReads != totalWrittenReads)
        {
            MD_LOGGER.warn("reads processed({}) vs written({}) mismatch", combinedStats.TotalReads, totalWrittenReads);
            chromosomeReaders.forEach(x -> x.recordWriter().logUnwrittenReads());
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

        if(mConfig.UnmapRegions.enabled())
        {
            MD_LOGGER.info("unmapped stats: {}", mConfig.UnmapRegions.stats().toString());
        }

        List<PerformanceCounter> combinedPerfCounters = chromosomeReaders.get(0).perfCounters();

        for(int i = 1; i < chromosomeReaders.size(); ++i)
        {
            List<PerformanceCounter> chrPerfCounters = chromosomeReaders.get(i).perfCounters();

            for(int j = 0; j < chrPerfCounters.size(); ++j)
            {
                combinedPerfCounters.get(j).merge(chrPerfCounters.get(j));
            }
        }

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

            List<Double> partitionLockTimes = Lists.newArrayList();
            partitionDataStore.partitions().forEach(x -> partitionLockTimes.add(x.totalLockTime()));
            Collections.sort(partitionLockTimes, Collections.reverseOrder());
            double nthTime = partitionLockTimes.size() >= 5 ? partitionLockTimes.get(4) : partitionLockTimes.get(partitionLockTimes.size() - 1);

            for(PartitionData partitionData : partitionDataStore.partitions())
            {
                double lockTime = partitionData.totalLockTime();

                if(lockTime > 0 && lockTime >= nthTime)
                {
                    MD_LOGGER.debug("partition({}) total lock-acquisition time({})",
                            partitionData.partitionStr(), format("%.3f", lockTime));
                }
            }
        }
        else
        {
            combinedPerfCounters.forEach(x -> x.logStats());
        }

        MD_LOGGER.info("Mark duplicates complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        MarkDuplicates markDuplicates = new MarkDuplicates(configBuilder);
        markDuplicates.run();
    }
}
