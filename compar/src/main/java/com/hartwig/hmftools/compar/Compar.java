package com.hartwig.hmftools.compar;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
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

    private BufferedWriter mDiffWriter;

    public Compar(final CommandLine cmd)
    {
        mConfig = new ComparConfig(cmd);

        mDiffWriter = null;
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            CMP_LOGGER.error("invalid config");
            return;
        }

        if(mConfig.SampleIds.isEmpty())
        {
            CMP_LOGGER.error("no samples specified");
            return;
        }

        if(mConfig.multiSample())
        {
            CMP_LOGGER.info("running comparison for {} sample(s)", mConfig.SampleIds.size());
        }

        initialiseOutputFiles();

        List<ComparTask> sampleTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            for(int i = 0; i < min(mConfig.SampleIds.size(), mConfig.Threads); ++i)
            {
                sampleTasks.add(new ComparTask(i, mConfig, mDiffWriter));
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
            ComparTask sampleTask = new ComparTask(0, mConfig, mDiffWriter);

            sampleTask.getSampleIds().addAll(mConfig.SampleIds);

            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        closeBufferedWriter(mDiffWriter);

        CMP_LOGGER.info("comparison complete");
    }

    private void initialiseOutputFiles()
    {
        try
        {
            String outputFile = mConfig.OutputDir;

            if(mConfig.singleSample())
            {
                outputFile += mConfig.SampleIds.get(0);
            }
            else
            {
                outputFile += "COMPAR_DIFFS";
            }

            if(mConfig.OutputId != null)
                outputFile += "." + mConfig.OutputId;

            outputFile += ".csv";

            mDiffWriter = createBufferedWriter(outputFile, false);

            if(mConfig.multiSample())
                mDiffWriter.write("SampleId,");

            mDiffWriter.write(Mismatch.header());

            mDiffWriter.newLine();
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to write compar output: {}", e.toString());
        }
    }

    public synchronized static void writeSampleMismatches(final BufferedWriter writer, final String sampleId, final List<Mismatch> mismatches)
    {
        if(mismatches.isEmpty() || writer == null)
            return;

        try
        {
            for(Mismatch mismatch : mismatches)
            {
                if(sampleId != null)
                    writer.write(String.format("%s,", sampleId));

                writer.write(mismatch.toCsv());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
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
