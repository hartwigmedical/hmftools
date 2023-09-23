package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public class RegionSlicer
{
    private final SliceConfig mConfig;

    public RegionSlicer(final ConfigBuilder configBuilder)
    {
        mConfig = new SliceConfig(configBuilder);
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

        for(ChrBaseRegion region : mConfig.SpecificChrRegions.Regions)
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
            String chromosome = entry.getKey();

            if(!HumanChromosome.contains(chromosome))
                continue;

            List<RemotePosition> remotePositions = entry.getValue();

            Collections.sort(remotePositions);

            remoteReadSlicers.add(new RemoteReadSlicer(entry.getKey(), remotePositions, mConfig, sliceWriter));
        }

        callableTasks = remoteReadSlicers.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
            System.exit(1);

        sliceWriter.close();

        BT_LOGGER.info("secondary slice complete");

        BT_LOGGER.info("Regions slice complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SliceConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        RegionSlicer regionSlicer = new RegionSlicer(configBuilder);
        regionSlicer.run();
    }
}
