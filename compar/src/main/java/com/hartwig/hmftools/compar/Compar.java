package com.hartwig.hmftools.compar;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.common.CommonUtils.initialiseFieldConfig;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FieldConfigFile;

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

        FieldConfig fieldConfig = initialiseFieldConfig(mConfig);
        fieldConfig.logProblems();
        if(fieldConfig.hasErrors())
        {
            System.exit(1);
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
                sampleTasks.add(new ComparTask(i, mConfig, fieldConfig, mWriter));
            }

            int taskIndex = 0;
            for(String sampleId : mConfig.SampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable<Void>> callableList = new ArrayList<>(sampleTasks);
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            ComparTask sampleTask = new ComparTask(0, mConfig, fieldConfig, mWriter);
            sampleTask.getSampleIds().addAll(mConfig.SampleIds);
            sampleTask.call();
        }

        mWriter.close();

        CMP_LOGGER.info("write field config file");
        try
        {
            Set<CategoryType> categories = buildComparers(mConfig).stream()
                    .map(c -> c.category())
                    .collect(Collectors.toSet());
            FieldConfigFile.write(FieldConfigFile.generateFileName(mConfig.OutputDir), fieldConfig, categories);
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("Could not write field config file", e);
            System.exit(1);
        }

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
