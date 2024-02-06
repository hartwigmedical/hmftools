package com.hartwig.hmftools.cup.liftover;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSomaticVcfFile;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class SnvLiftover
{
    private final LiftoverConfig mConfig;
    private final List<String> mSampleIds;
    private final GenomeLiftoverCache mMappingCache;
    private final int mThreads;

    @Deprecated
    public static final String LIFTOVER_FILE = ".snv_liftover.csv";

    public SnvLiftover(final ConfigBuilder configBuilder)
    {
        mConfig = new LiftoverConfig(configBuilder);

        if(configBuilder.hasValue(SAMPLE))
        {
            mSampleIds = Lists.newArrayList(configBuilder.getValue(SAMPLE));
        }
        else
        {
            mSampleIds = loadSampleIdsFile(configBuilder);
        }
        mThreads = parseThreads(configBuilder);

        mMappingCache = new GenomeLiftoverCache(true);
    }

    public void run()
    {
        if(mConfig.OutputDir == null || mConfig.SampleVcfDir == null || mSampleIds.isEmpty())
        {
            CUP_LOGGER.error("invalid config");
            System.exit(1);
        }

        CUP_LOGGER.info("converting positions for {} samples", mSampleIds.size());

        List<VcfPositionConverter> sampleTasks = Lists.newArrayList();

        for(String sampleId : mSampleIds)
        {
            String purpleDir = mConfig.SampleVcfDir.replaceAll("\\*", sampleId);
            String vcfFile = purpleSomaticVcfFile(purpleDir, sampleId);
            VcfPositionConverter vcfTask = new VcfPositionConverter(sampleId, vcfFile, mMappingCache, mConfig);
            sampleTasks.add(vcfTask);
        }

        if(mThreads > 1)
        {
            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mThreads);
        }
        else
        {
            sampleTasks.forEach(x -> x.call());
        }

        CUP_LOGGER.info("conversion complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        LiftoverConfig.addOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        SnvLiftover snvLiftover = new SnvLiftover(configBuilder);
        snvLiftover.run();
    }
}
