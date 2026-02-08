package com.hartwig.hmftools.isofox.loader;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.loader.DataLoadType.GENE_EXPRESSION;
import static com.hartwig.hmftools.isofox.loader.DataLoadType.NOVEL_JUNCTION;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.expression.cohort.CohortGenePercentiles;
import com.hartwig.hmftools.isofox.novel.AltSjCohortCache;

import org.jetbrains.annotations.NotNull;

public class IsofoxDataLoader
{
    private final DataLoaderConfig mConfig;

    public IsofoxDataLoader(final ConfigBuilder configBuilder)
    {
        mConfig = new DataLoaderConfig(configBuilder);
    }

    public boolean load()
    {
        if(mConfig.SampleIds.isEmpty())
        {
            ISF_LOGGER.error("no sample IDs configured");
            return false;
        }

        List<SampleLoaderTask> sampleTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            for(int i = 0; i < min(mConfig.SampleIds.size(), mConfig.Threads); ++i)
            {
                sampleTasks.add(new SampleLoaderTask(i, mConfig));
            }

            int taskIndex = 0;
            for(String sampleId : mConfig.SampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            SampleLoaderTask sampleTask = new SampleLoaderTask(0, mConfig);

            sampleTask.getSampleIds().addAll(mConfig.SampleIds);

            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        ISF_LOGGER.info("Isofox data loading complete");

        return true;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        DataLoaderConfig.registerConfig(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        IsofoxDataLoader dataLoader = new IsofoxDataLoader(configBuilder);

        if(!dataLoader.load())
        {
            ISF_LOGGER.info("Isofox data loading failed");
            System.exit(1);
        }

    }
}
