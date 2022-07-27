package com.hartwig.hmftools.svprep.tools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class HighDepthFinder
{
    private final HighDepthConfig mConfig;
    private final BufferedWriter mWriter;

    public HighDepthFinder(final CommandLine cmd)
    {
        mConfig = new HighDepthConfig(cmd);
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
            if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chromosome.toString())))
                continue;

            HighDepthTask depthTask = new HighDepthTask(chromosome.toString(), mConfig, mWriter);
            depthTasks.add(depthTask);
        }

        final List<Callable> callableList = depthTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("high depth discovery complete");

        // write output VCF
    }

    private BufferedWriter initialiseWriter(final String filename)
    {
        SV_LOGGER.info("writing output to {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome,PosStart,PosEnd,BaseDepth");
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
                writer.write(format("%s,%d,%d,%d", region.Region.Chromosome, region.Region.start(), region.Region.end(), region.Depth));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write region: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        HighDepthConfig.addOptions(options);
        addSpecificChromosomesRegionsConfig(options);
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        HighDepthFinder highDepthFinder = new HighDepthFinder(cmd);
        highDepthFinder.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
