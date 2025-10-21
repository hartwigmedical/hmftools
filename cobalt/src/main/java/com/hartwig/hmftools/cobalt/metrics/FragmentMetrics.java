package com.hartwig.hmftools.cobalt.metrics;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.APP_NAME;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FragmentMetrics
{
    private final MetricsConfig mConfig;

    public FragmentMetrics(final ConfigBuilder configBuilder)
    {
        mConfig = new MetricsConfig(configBuilder);
    }

    public void run() throws IOException
    {
        CB_LOGGER.info("sample({}) starting fragment GC metrics", mConfig.SampleId);
        long startTimeMs = System.currentTimeMillis();
        ListMultimap<Chromosome, Partition> partitionsByChromosome = mConfig.createPartitions();
        Queue<Partition> partitions = new ConcurrentLinkedQueue<>(partitionsByChromosome.values());
        CB_LOGGER.info("splitting {} partitions across {} threads", partitions.size(), mConfig.Threads);

        List<Thread> workers = new ArrayList<>();
        for(int i = 0; i < mConfig.Threads; ++i)
        {
            PartitionReader partitionReader = new PartitionReader(mConfig, partitions);
            partitionReader.start();
            workers.add(partitionReader);
        }
        if(!runThreadTasks(workers))
        {
            System.exit(1);
        }
        CB_LOGGER.debug("all partitions complete");

        List<WindowStatistics> results = new ArrayList<>(partitionsByChromosome.size());
        partitionsByChromosome.keySet().forEach(chromosome ->
        {
            List<Partition> partitionsForChromosome = partitionsByChromosome.get(chromosome);
            partitionsForChromosome.forEach(partition -> partition.TargetRegions.forEach(region -> results.add(region.statistics())));
        });
        String filename = WindowStatisticsFile.fileName(mConfig.OutputDir, mConfig.SampleId);
        CB_LOGGER.info("writing fragment statistics to: {}", filename);
        WindowStatisticsFile.write(filename, results);
        CB_LOGGER.info("FragmentGcMetrics complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(final String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        MetricsConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        FragmentMetrics fragmentMetrics = new FragmentMetrics(configBuilder);
        fragmentMetrics.run();
    }
}
