package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.NotNull;

// TODO(m_cooper): BED and blacklist region repeater mask mismatch. Zeroes in the BED file.
public class Annotate
{
    private final AnnotateConfig mConfig;

    public Annotate(final ConfigBuilder configBuilder)
    {
        mConfig = new AnnotateConfig(configBuilder);
    }

    public void run()
    {
        if(mConfig.BamFile == null)
        {
            MD_LOGGER.error("no BAM file specified");
            System.exit(1);
        }

        if(mConfig.BlacklistRegionRepeatMaskerFile == null)
        {
            MD_LOGGER.error("no blacklist region repeat masker CSV file specified");
            System.exit(1);
        }

        final List<AnnotatedBlacklistRegion> blacklistRegions =
                BlacklistRepeatMaskerReader.readFromFile(mConfig.BlacklistRegionRepeatMaskerFile);

        final List<IJobRegion> jobRegions = new ArrayList<>();
        for(AnnotatedBlacklistRegion blacklistRegion : blacklistRegions)
        {
            jobRegions.addAll(blacklistRegion.getJobRegions());
        }

        // TODO(m_cooper): Overlapping
        // TODO(m_cooper): 9
        // TODO(m_cooper): [68999357, 69000000]
        // TODO(m_cooper): [69000000, 69000299]
        // TODO(m_cooper): PosEnd included?

        // TODO(m_cooper): Read entirely contained?

        final ArrayBlockingQueue<IJobRegion> jobs = new ArrayBlockingQueue<>(jobRegions.size(), true, jobRegions);
        final List<AnnotateConsumer> annotateConsumers = new ArrayList<>();
        for(int i = 0; i < Math.max(mConfig.Threads, 1); i++)
        {
            final AnnotateConsumer annotateConsumer = new AnnotateConsumer(mConfig, jobs);
            annotateConsumers.add(annotateConsumer);
        }

        final List<Callable> callableList = annotateConsumers.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        writeAnnotatedBedFile(blacklistRegions);

        MD_LOGGER.info("annotate complete");
    }

    private void writeAnnotatedBedFile(@NotNull final List<AnnotatedBlacklistRegion> blacklistRegions)
    {
        MD_LOGGER.info("Writing output to {}.", mConfig.OutputFile);
        try
        {
            BufferedWriter writer = new BufferedWriter(new FileWriter(mConfig.OutputFile));
            writer.write(AnnotatedBlacklistRegion.CSV_HEADER);
            writer.newLine();
            for(AnnotatedBlacklistRegion blacklistRegion : blacklistRegions)
            {
                for(String line : blacklistRegion.getCSVLines())
                {
                    writer.write(line);
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(Exception e)
        {
            MD_LOGGER.error("An exception was raised while writing the output to {}: {}", mConfig.OutputFile, e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        AnnotateConfig.addConfig(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);
        // TODO(m_cooper): Fill in.
        // logVersion();

        Annotate annotate = new Annotate(configBuilder);
        annotate.run();
    }
}
