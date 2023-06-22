package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.PartitionTask.partitionChromosome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

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

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            allRegions.addAll(partitionChromosome(chromosomeStr, mConfig.RefGenVersion, mConfig.SpecificRegions, mConfig.PartitionSize));
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

        BT_LOGGER.info("all regions complete, totalReads({})", combinedStats.totalReads());

        combinedStats.coverageMetrics().finalise(mConfig.ExcludeZeroCoverage);
        MetricsWriter.writeResults(combinedStats, mConfig);
        BT_LOGGER.info("final stats: {}", combinedStats.coverageMetrics());

        BT_LOGGER.info("all regions complete, totalReads({}) stats: {}", combinedStats.totalReads(), combinedStats.coverageMetrics());

        if(mConfig.PerfDebug)
        {
            combinedStats.perfCounter().logIntervalStats(10);
            combinedStats.perfCounter().logStats();
        }

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        BT_LOGGER.info("BamMetrics complete, mins({})", format("%.3f", timeTakeMins));
    }

    public static void main(@NotNull final String[] args)
    {
        final VersionInfo version = new VersionInfo("bam-tools.version");
        BT_LOGGER.info("BamTools version: {}", version.version());

        ConfigBuilder configBuilder = new ConfigBuilder();
        MetricsConfig.addConfig(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logItems();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        BamMetrics bamMetrtics = new BamMetrics(configBuilder);
        bamMetrtics.run();
    }
}
