package com.hartwig.hmftools.cup.liftover;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.LIFTOVER_MAPPING_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSomaticVcfFile;
import static com.hartwig.hmftools.cup.liftover.LiftoverConfig.SAMPLE;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class SnvLiftover
{
    private final LiftoverConfig mConfig;
    private final List<String> mSampleIds;
    private final GenomeLiftoverCache mMappingCache;
    private final int mThreads;

    @Deprecated
    public static final String LIFTOVER_FILE = ".snv_liftover.csv";

    public SnvLiftover(final CommandLine cmd)
    {
        mConfig = new LiftoverConfig(cmd);

        if(cmd.hasOption(SAMPLE))
        {
            mSampleIds = Lists.newArrayList(cmd.getOptionValue(SAMPLE));
        }
        else
        {
            mSampleIds = loadSampleIdsFile(cmd);
        }
        mThreads = parseThreads(cmd);

        mMappingCache = new GenomeLiftoverCache();

        if(cmd.hasOption(LIFTOVER_MAPPING_FILE))
        {
            if(!mMappingCache.loadFile(cmd.getOptionValue(LIFTOVER_MAPPING_FILE)))
                System.exit(1);
        }
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        LiftoverConfig.addOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        SnvLiftover snvLiftover = new SnvLiftover(cmd);
        snvLiftover.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
