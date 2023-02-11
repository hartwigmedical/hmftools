package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class RegionSlicer
{
    private final SliceConfig mConfig;

    public RegionSlicer(final CommandLine cmd)
    {
        mConfig = new SliceConfig(cmd);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        BT_LOGGER.info("sample({}) starting BamSlicer", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        SliceWriter sliceWriter = new SliceWriter(mConfig);
        ReadCache readCache = new ReadCache(mConfig);

        List<RegionBamSlicer> regionBamSlicers = Lists.newArrayList();

        for(ChrBaseRegion region : mConfig.SpecificRegions)
        {
            regionBamSlicers.add(new RegionBamSlicer(region, mConfig, readCache, sliceWriter));
        }

        BT_LOGGER.info("splitting {} regions across {} threads", regionBamSlicers.size(), mConfig.Threads);

        List<Callable> callableTasks = regionBamSlicers.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
            System.exit(1);

        BT_LOGGER.info("initial slice complete");

        List<RemoteReadSlicer> remoteReadSlicers = Lists.newArrayList();

        for(Map.Entry<String,List<RemotePosition>> entry : readCache.chrRemotePositions().entrySet())
        {
            remoteReadSlicers.add(new RemoteReadSlicer(entry.getKey(), entry.getValue(), mConfig, sliceWriter));

        }

        callableTasks = remoteReadSlicers.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
            System.exit(1);

        sliceWriter.close();

        BT_LOGGER.info("secondary slice complete");

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        BT_LOGGER.info("Regions slice complete, mins({})", format("%.3f", timeTakeMins));
    }

    public static void main(@NotNull final String[] args)
    {
        final VersionInfo version = new VersionInfo("bam-tools.version");
        BT_LOGGER.info("BamTools version: {}", version.version());

        final Options options = SliceConfig.createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            RegionSlicer regionSlicer = new RegionSlicer(cmd);
            regionSlicer.run();
        }
        catch(ParseException e)
        {
            BT_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("RegionSlicer", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
