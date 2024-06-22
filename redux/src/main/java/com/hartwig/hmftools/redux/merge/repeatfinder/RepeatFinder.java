package com.hartwig.hmftools.redux.merge.repeatfinder;

import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.merge.repeatfinder.RepeatFinderConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;

public class RepeatFinder
{
    private static final String OUTPUT_TSV_HEADER = "Chromosome\tPosStart\tPosEnd\tRepeatBases";
    private static final String BED_ANNOTATION_TSV_HEADER = "WithinBedRegions\tPartialOverlapsBedRegions";

    private final RepeatFinderConfig mConfig;

    public RepeatFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new RepeatFinderConfig(configBuilder);
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        RepeatFinderConfig.addConfig(configBuilder);

        addLoggingOptions(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        RepeatFinder repeatFinder = new RepeatFinder(configBuilder);
        repeatFinder.run();
    }

    public static synchronized void writeRepeat(final BufferedWriter writer, final RepeatInfo repeatInfo)
    {
        try
        {
            writer.write(repeatInfo.getTSVFragment());
            writer.newLine();
        }
        catch(IOException e)
        {
            MD_LOGGER.error("Failed to write repeat to output file: {}", e.toString());
            System.exit(1);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            MD_LOGGER.info("Writing output to {}", mConfig.OutputFile);
            BufferedWriter writer = createBufferedWriter(mConfig.OutputFile, false);

            writer.write(OUTPUT_TSV_HEADER);
            if(mConfig.BedFile != null)
            {
                writer.write("\t" + BED_ANNOTATION_TSV_HEADER);
            }
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            MD_LOGGER.error("Failed to initialise output file: {}", e.toString());
            System.exit(1);
            return null;
        }
    }

    public void run()
    {
        // See docstring of `RepeatFinderTask.processChromosome` for an explanation of this requirement.
        if(mConfig.MinRepeatBases < 4)
        {
            MD_LOGGER.error("Value of -min_repeat_bases must be at least four.");
            System.exit(1);
        }

        Map<String, List<BaseRegion>> chrRegionsMap = null;
        if(mConfig.BedFile != null)
        {
            chrRegionsMap = loadChrBaseRegions(mConfig.BedFile);
            chrRegionsMap.values().stream().forEach(Collections::sort);
        }

        BufferedWriter writer = initialiseWriter();

        List<RepeatFinderTask> finderTasks = Lists.newArrayList();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());
            RepeatFinderTask finderTask;
            if(chrRegionsMap == null)
            {
                finderTask = new RepeatFinderTask(mConfig, chrStr, writer);
                finderTasks.add(finderTask);
                continue;
            }

            List<BaseRegion> chrRegions = chrRegionsMap.get(chrStr);
            finderTask = new RepeatFinderTask(mConfig, chrStr, writer, chrRegions);
            finderTasks.add(finderTask);
        }

        List<Callable> callableList = finderTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        closeBufferedWriter(writer);

        MD_LOGGER.info("repeat finder complete");
    }
}
