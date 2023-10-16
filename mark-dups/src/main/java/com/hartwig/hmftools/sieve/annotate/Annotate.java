package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class Annotate
{
    private final AnnotateConfig mConfig;

    public Annotate(final ConfigBuilder configBuilder)
    {
        mConfig = new AnnotateConfig(configBuilder);
    }

    private static String regionToTSVFragement(final ChrBaseRegion region)
    {
        return region.chromosome() + '\t' + region.start() + '\t' + region.end();
    }

    public static synchronized void writeRecord(
            final BufferedWriter outputWriter,
            final ChrBaseRegion region,
            final RefGenomeRegionAnnotations refAnnotations,
            final HighDepthCounts counts)
    {
        try
        {
            outputWriter.write(regionToTSVFragement(region) + '\t' + refAnnotations.getTSVFragment() + '\t' + counts.getTSVFragment());
            outputWriter.newLine();
        }
        catch(IOException e)
        {
            MD_LOGGER.error("An exception was raised while writing a record to the output file: {}", e.toString());
            System.exit(1);
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        AnnotateConfig.addConfig(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);
        // TODO(m_cooper): Fill in.
        // logVersion();

        Annotate annotate = new Annotate(configBuilder);
        annotate.run();
    }

    public void run()
    {
        if(mConfig.BamFile == null)
        {
            MD_LOGGER.error("no BAM file specified");
            System.exit(1);
        }

        if(mConfig.HighDepthFile == null)
        {
            MD_LOGGER.error("no high depth file specified");
            System.exit(1);
        }

        BufferedWriter outputWriter = createOutputWriter();
        final List<ChrBaseRegion> regions = HighDepthReader.readFromFile(mConfig.HighDepthFile);

        final ArrayBlockingQueue<ChrBaseRegion> jobs = new ArrayBlockingQueue<>(regions.size(), true, regions);
        final List<HighDepthCountConsumer> annotateConsumers = new ArrayList<>();
        for(int i = 0; i < Math.max(mConfig.Threads, 1); i++)
        {
            final HighDepthCountConsumer consumer = new HighDepthCountConsumer(mConfig, jobs, outputWriter);
            annotateConsumers.add(consumer);
        }

        final List<Callable> callableList = annotateConsumers.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        try
        {
            outputWriter.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error("An exception was raised while closing the output file {}: {}", mConfig.OutputFile, e.toString());
            System.exit(1);
        }

        MD_LOGGER.info("annotate complete");
    }

    private BufferedWriter createOutputWriter()
    {
        MD_LOGGER.info("Writing output to {}.", mConfig.OutputFile);
        try
        {
            BufferedWriter writer = new BufferedWriter(new FileWriter(mConfig.OutputFile));
            writer.write("Chromosome\tPosStart\tPosEnd\t" + RefGenomeRegionAnnotations.TSV_HEADER + '\t' + HighDepthCounts.TSV_HEADER);
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            MD_LOGGER.error("An exception was raised while initialising a BufferedWriter for the output to {}: {}", mConfig.OutputFile, e.toString());
            System.exit(1);
        }

        return null;
    }
}
