package com.hartwig.hmftools.compar;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class Compar
{
    private final ComparConfig mConfig;
    private final MismatchWriter mWriter;

    public Compar(final CommandLine cmd)
    {
        mConfig = new ComparConfig(cmd);
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

        List<ComparTask> sampleTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
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

            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            ComparTask sampleTask = new ComparTask(0, mConfig, mWriter);

            sampleTask.getSampleIds().addAll(mConfig.SampleIds);

            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        mWriter.close();

        CMP_LOGGER.info("comparison complete");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        ComparConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        Compar compar = new Compar(cmd);
        compar.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
