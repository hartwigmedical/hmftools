package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BamMetrics
{
    private final MetricsConfig mConfig;

    public BamMetrics(final ConfigBuilder configBuilder)
    {
        mConfig = new MetricsConfig(configBuilder);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        BT_LOGGER.info("sample({}) starting bam metrics", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        List<ChrBaseRegion> allRegions = Lists.newArrayList();

        if(!mConfig.TargetRegions.isEmpty())
        {
            int targetRegionCount = mConfig.TargetRegions.values().stream().mapToInt(x -> x.size()).sum();

            int targetRegionBaseCount = mConfig.TargetRegions.values().stream()
                    .mapToInt(x -> x.stream().mapToInt(y -> y.baseLength()).sum()).sum();

            BT_LOGGER.info("capturing data for {} target regions, total base count({})",
                    targetRegionCount, targetRegionBaseCount);
        }

        if(!mConfig.OnlyTargetRegions)
        {
            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

                if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                    continue;

                allRegions.addAll(partitionChromosome(
                        chromosomeStr, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize));
            }
        }
        else
        {
            for(Map.Entry<String,List<BaseRegion>> entry : mConfig.TargetRegions.entrySet())
            {
                entry.getValue().forEach(x -> allRegions.add(new ChrBaseRegion(entry.getKey(), x.start(), x.end())));
            }
        }

        BT_LOGGER.info("splitting {} regions across {} threads", allRegions.size(), mConfig.Threads);

        Queue<PartitionTask> partitions = new ConcurrentLinkedQueue<>();

        int taskId = 0;
        for(int i = 0; i < allRegions.size(); ++i)
        {
            partitions.add(new PartitionTask(allRegions.get(i), taskId++));
        }

        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(allRegions.size(), mConfig.Threads); ++i)
        {
            workers.add(new PartitionThread(mConfig, partitions, combinedStats));
        }

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                BT_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        BT_LOGGER.info("all regions complete");

        combinedStats.coverageMetrics().finalise(mConfig.ExcludeZeroCoverage);
        MetricsWriter.writeResults(combinedStats, mConfig);

        BT_LOGGER.info("totalReads({}) stats: {}", combinedStats.readCounts().TotalReads, combinedStats.coverageMetrics());

        if(mConfig.PerfDebug)
        {
            combinedStats.perfCounter().logIntervalStats(10);
            combinedStats.perfCounter().logStats();
        }

        BT_LOGGER.info("BamMetrics complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        MetricsConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BamMetrics bamMetrtics = new BamMetrics(configBuilder);
        bamMetrtics.run();
    }
}
