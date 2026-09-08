package com.hartwig.hmftools.bamtools.depth;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CRAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.NotNull;

public class HighDepthFinder
{
    private final FinderConfig mConfig;
    private final BufferedWriter mWriter;

    public HighDepthFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new FinderConfig(configBuilder);

        String outputFile;

        if(mConfig.OutputFile != null)
        {
            outputFile = mConfig.OutputFile;
        }
        else
        {
            int bamCramIndex = mConfig.BamFile.endsWith(CRAM_EXTENSION) ?
                    mConfig.BamFile.indexOf(CRAM_EXTENSION) : mConfig.BamFile.indexOf(BAM_EXTENSION);

            outputFile = mConfig.BamFile.substring(0, bamCramIndex) + "." + HIGH_DEPTH_FILE_ID + TSV_EXTENSION;
        }

        mWriter = initialiseWriter(outputFile);
    }

    public void run()
    {
        if(mConfig.BamFile == null)
        {
            BT_LOGGER.error("no BAM file specified");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        List<HighDepthTask> depthTasks = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chrStr)))
                continue;

            HighDepthTask depthTask = new HighDepthTask(chrStr, mConfig, mWriter);
            depthTasks.add(depthTask);
        }

        List<Callable<Void>> callableList = depthTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        closeBufferedWriter(mWriter);

        BT_LOGGER.info("High depth finder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    protected static final String FLD_BASE_DEPTH_MAX = "BaseDepthMax";
    protected static final String FLD_BASE_DEPTH_MIN = "BaseDepthMin";
    protected static final String FLD_BASE_DEPTH_AVG = "BaseDepthAvg";
    protected static final String HIGH_DEPTH_FILE_ID = "high_depth";

    private BufferedWriter initialiseWriter(final String filename)
    {
        BT_LOGGER.info("writing output to {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME).add(FLD_POS_START).add(FLD_POS_END);
            sj.add(FLD_BASE_DEPTH_MIN).add(FLD_BASE_DEPTH_MAX).add(FLD_BASE_DEPTH_AVG);
            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to initialise writer: {}", e.toString());
        }

        return null;
    }

    public synchronized static void writeHighDepthRegions(final BufferedWriter writer, final List<HighDepthRegion> regions)
    {
        if(writer == null)
            return;

        try
        {
            for(HighDepthRegion region : regions)
            {
                if(region.DepthAvg < region.DepthMin)
                    continue;

                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(region.Chromosome);
                sj.add(String.valueOf(region.start()));
                sj.add(String.valueOf(region.end()));
                sj.add(String.valueOf(region.DepthMin));
                sj.add(String.valueOf(region.DepthMax));
                sj.add(String.valueOf(region.DepthAvg));
                writer.write(sj.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to write region: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        FinderConfig.addConfig(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        HighDepthFinder highDepthFinder = new HighDepthFinder(configBuilder);
        highDepthFinder.run();
    }
}
