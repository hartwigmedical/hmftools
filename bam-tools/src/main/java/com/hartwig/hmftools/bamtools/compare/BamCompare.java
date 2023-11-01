package com.hartwig.hmftools.bamtools.compare;

import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public class BamCompare
{
    private final CompareConfig mConfig;

    public BamCompare(final ConfigBuilder configBuilder)
    {
        mConfig = new CompareConfig(configBuilder);
    }

    public void run()
    {
        BT_LOGGER.info("starting bam comparison", mConfig.OutputFile);

        long startTimeMs = System.currentTimeMillis();

        List<ChrBaseRegion> allRegions = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            List<ChrBaseRegion> partitions = partitionChromosome(
                    chromosomeStr, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize);

            allRegions.addAll(partitions);
        }

        ReadWriter readWriter = new ReadWriter(mConfig);

        if(!readWriter.initialised())
            System.exit(1);

        BT_LOGGER.info("splitting {} regions across {} threads", allRegions.size(), mConfig.Threads);

        Queue<PartitionTask> partitions = new ConcurrentLinkedQueue<>();

        int taskId = 0;
        for(int i = 0; i < allRegions.size(); ++i)
        {
            partitions.add(new PartitionTask(allRegions.get(i), taskId++));
        }

        List<PartitionThread> partitionTasks = Lists.newArrayList();
        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(allRegions.size(), mConfig.Threads); ++i)
        {
            PartitionThread partitionThread = new PartitionThread(mConfig, partitions, readWriter);
            partitionTasks.add(partitionThread);
            workers.add(partitionThread);
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

        Statistics combinedStats = new Statistics();
        partitionTasks.forEach(x -> combinedStats.merge(x.stats()));

        readWriter.close();

        BT_LOGGER.info("summary: reads(ref={} new={}) diffs({})",
                combinedStats.RefReadCount, combinedStats.NewReadCount, combinedStats.DiffCount);

        BT_LOGGER.info("BamCompare complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BamCompare bamCompare = new BamCompare(configBuilder);
        bamCompare.run();
    }
}
