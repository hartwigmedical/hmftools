package com.hartwig.hmftools.cup.liftover;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSomaticVcfFile;

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

public class SnvLiftover
{
    private final String mOutputDir;
    private final String mSampleVcfDir;
    private final List<String> mSampleIds;
    private final CoordMappingCache mMappingCache;
    private final int mThreads;

    // config
    private static final String SAMPLE = "sample";
    private static final String SAMPLE_VCF_DIR = "sample_vcf_dir";
    private static final String MAPPING_FILE = "mapping_file";

    public SnvLiftover(final CommandLine cmd)
    {
        if(cmd.hasOption(SAMPLE))
        {
            mSampleIds = Lists.newArrayList(cmd.getOptionValue(SAMPLE));
        }
        else
        {
            mSampleIds = loadSampleIdsFile(cmd);
        }
        mOutputDir = parseOutputDir(cmd);
        mSampleVcfDir = cmd.getOptionValue(SAMPLE_VCF_DIR);
        mThreads = parseThreads(cmd);

        mMappingCache = new CoordMappingCache();

        try
        {
            mMappingCache.loadFile(cmd.getOptionValue(MAPPING_FILE));
        }
        catch(Exception e)
        {
            System.exit(1);
        }
    }

    public void run()
    {
        if(mOutputDir == null || mSampleVcfDir == null || mSampleIds.isEmpty())
        {
            CUP_LOGGER.error("invalid config");
            System.exit(1);
        }

        CUP_LOGGER.info("converting positions for {} samples", mSampleIds.size());

        List<VcfPositionConverter> sampleTasks = Lists.newArrayList();

        for(String sampleId : mSampleIds)
        {
            String purpleDir = mSampleVcfDir.replaceAll("\\*", sampleId);
            String vcfFile = purpleSomaticVcfFile(purpleDir, sampleId);
            VcfPositionConverter vcfTask = new VcfPositionConverter(sampleId, vcfFile, mOutputDir, mMappingCache);
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

        addSampleIdFile(options);
        addOutputOptions(options);
        addThreadOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Sample ID");
        options.addOption(SAMPLE_VCF_DIR, true, "Path to sample VCF(s)");
        options.addOption(MAPPING_FILE, true, "Coordinate mapping file");

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
