package com.hartwig.hmftools.compar;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class Compar
{
    private final ComparConfig mConfig;
    private final MismatchWriter mWriter;

    public Compar(final ConfigBuilder configBuilder)
    {
        mConfig = new ComparConfig(configBuilder);
        mWriter = new MismatchWriter(mConfig);
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            CMP_LOGGER.error("invalid config");
            System.exit(1);
        }

        if(mConfig.SampleIds.isEmpty())
        {
            CMP_LOGGER.error("no samples specified");
            System.exit(1);
        }

        if(mConfig.multiSample())
        {
            CMP_LOGGER.info("running comparison for {} sample(s)", mConfig.SampleIds.size());
        }

        if(!mWriter.initialiseOutputFiles())
        {
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        if(mConfig.Threads > 1)
        {
            List<ComparTask> sampleTasks = Lists.newArrayList();

            for(int i = 0; i < min(mConfig.SampleIds.size(), mConfig.Threads); ++i)
            {
                sampleTasks.add(new ComparTask(i, mConfig, mWriter));
            }

            int taskIndex = 0;
            for(String sampleId : mConfig.SampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable> callableList = new ArrayList<>(sampleTasks);
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            ComparTask sampleTask = new ComparTask(0, mConfig, mWriter);
            sampleTask.getSampleIds().addAll(mConfig.SampleIds);
            sampleTask.call();
        }

        mWriter.close();

        if(mConfig.multiSample())
        {
            CMP_LOGGER.info("comparison of {} samples complete, mins({})", mConfig.SampleIds.size(), runTimeMinsStr(startTimeMs));
        }
        else
        {
            CMP_LOGGER.info("comparison complete");
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Compar");
        ComparConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        Compar compar = new Compar(configBuilder);
        compar.run();
    }
}
