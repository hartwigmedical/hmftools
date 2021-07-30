package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;
import static com.hartwig.hmftools.telo.TeloUtils.createPartitions;

import java.io.IOException;
import java.time.Instant;
import java.time.Duration;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;

public class TeloApplication
{
    private final TeloConfig mConfig;

    private final BamRecordWriter mWriter;

    public TeloApplication(final Options options, final String... args) throws ParseException
    {
        VersionInfo versionInfo = new VersionInfo("telo.version");
        TE_LOGGER.info("Telo version: {}", versionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        mConfig = new TeloConfig(cmd);

        mWriter = new BamRecordWriter(mConfig);
    }

    private void run()
    {
        if(!mConfig.isValid())
        {
            TE_LOGGER.error(" invalid config, exiting");
            System.exit(1);
        }

        TE_LOGGER.info("starting telomeric analysis");
        Instant start = Instant.now();

        processBam();

        mWriter.close();

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        TE_LOGGER.info("Telo run complete, time taken: {}m {}s",  seconds / 60, seconds % 60);
    }

    private void processBam()
    {
        final ConcurrentMap<Thread, SamReader> threadSamReaders = new ConcurrentHashMap<>();
        final ConcurrentMap<String, ReadGroup> incompleteReadGroups = new ConcurrentHashMap<>();;
        final List<BamReader> bamReaders = Lists.newArrayList();
        final List<Callable> callableList = Lists.newArrayList();

        Arrays.stream(HumanChromosome.values())
                .map(chromosome -> chromosome.toString());

        List<BaseRegion> partitions = createPartitions(mConfig);

        for(BaseRegion partition : partitions)
        {
            BamReader bamReader = new BamReader(mConfig, mWriter, threadSamReaders, incompleteReadGroups);
            bamReader.setBaseRegion(partition);

            bamReaders.add(bamReader);
            callableList.add(bamReader);
        }

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
        {
            TE_LOGGER.error("BAM file processing failed");
            return;
        }

        TE_LOGGER.info("initial BAM file processing complete");

        BamReader finalBamReader = new BamReader(mConfig, mWriter, threadSamReaders, incompleteReadGroups);
        finalBamReader.findTelomereContent();

        TE_LOGGER.info("writing final {} incomplete groups", incompleteReadGroups.size());

        incompleteReadGroups.values().forEach(x -> mWriter.writeReadGroup(x));
    }

    public static void main(final String... args) throws IOException
    {
        final Options options = TeloConfig.createOptions();

        try
        {
            TeloApplication application = new TeloApplication(options, args);
            application.run();
        }
        catch(ParseException e)
        {
            TE_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CountBamLinesApplication", options);
            System.exit(1);
        }
    }

    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
