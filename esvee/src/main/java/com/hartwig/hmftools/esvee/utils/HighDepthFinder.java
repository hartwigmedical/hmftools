package com.hartwig.hmftools.esvee.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.NotNull;

public class HighDepthFinder
{
    private final HighDepthConfig mConfig;
    private final BufferedWriter mWriter;

    public HighDepthFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new HighDepthConfig(configBuilder);
        mWriter = initialiseWriter(mConfig.OutputFile);
    }

    public void run()
    {
        if(mConfig.BamFile == null)
        {
            SV_LOGGER.error("no BAM file specified");
            System.exit(1);
        }

        List<HighDepthTask> depthTasks = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chrStr)))
                continue;

            HighDepthTask depthTask = new HighDepthTask(chrStr, mConfig, mWriter);
            depthTasks.add(depthTask);
        }

        final List<Callable> callableList = depthTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("high depth discovery complete");
    }

    private BufferedWriter initialiseWriter(final String filename)
    {
        SV_LOGGER.info("writing output to {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome\tPosStart\tPosEnd\tBaseDepthMin\tBaseDepthMax");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to initialise writer: {}", e.toString());
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
                writer.write(format("%s\t%d\t%d\t%d\t%d",
                        region.Chromosome, region.start(), region.end(), region.DepthMin, region.DepthMax));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write region: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        HighDepthConfig.addConfig(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        HighDepthFinder highDepthFinder = new HighDepthFinder(configBuilder);
        highDepthFinder.run();
    }
}
