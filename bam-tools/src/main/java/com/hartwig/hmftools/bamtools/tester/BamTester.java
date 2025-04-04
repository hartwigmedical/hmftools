package com.hartwig.hmftools.bamtools.tester;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;

import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.checker.CheckConfig;
import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskQueue;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BamTester
{
    private final TesterConfig mConfig;

    public BamTester(final ConfigBuilder configBuilder)
    {
        mConfig = new TesterConfig(configBuilder);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        BT_LOGGER.info("running BAM tester");

        List<ChrBaseRegion> partitionRegions = PartitionTask.splitRegionsIntoPartitions(
                mConfig.BamFile, mConfig.RefGenomeFile, mConfig.Threads, mConfig.SpecificChrRegions, mConfig.PartitionSize);

        List<PartitionTester> partitionTesters = Lists.newArrayList();
        List<Thread> threadTasks = Lists.newArrayList();

        Queue<ChrBaseRegion> regionsQueue = new ConcurrentLinkedQueue<>();
        regionsQueue.addAll(partitionRegions);

        TaskQueue taskQueue = new TaskQueue(regionsQueue, "regions", 100);

        for(int i = 0; i < mConfig.Threads; ++i)
        {
            PartitionTester partitionTester = new PartitionTester(mConfig, taskQueue);
            partitionTesters.add(partitionTester);
            threadTasks.add(partitionTester);
        }

        BT_LOGGER.debug("splitting {} regions across {} threads", partitionRegions.size(), mConfig.Threads);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        long totalReads = partitionTesters.stream().mapToLong(x -> x.totalReads()).sum();

        BT_LOGGER.info("processed {} reads", totalReads);

        logPerformanceStats(partitionTesters);

        BT_LOGGER.info("BAM tester complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void logPerformanceStats(final List<PartitionTester> partitionTesters)
    {
        List<PerformanceCounter> combinedPerfCounters = partitionTesters.get(0).perfCounters();

        for(int i = 1; i < partitionTesters.size(); ++i)
        {
            List<PerformanceCounter> threadPerfCounters = partitionTesters.get(i).perfCounters();

            for(int j = 0; j < combinedPerfCounters.size(); ++j)
            {
                combinedPerfCounters.get(j).merge(threadPerfCounters.get(j));
            }
        }

        for(PerformanceCounter perfCounter : combinedPerfCounters)
        {
            perfCounter.logStats();
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        TesterConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BamTester bamTester = new BamTester(configBuilder);
        bamTester.run();
    }
}
