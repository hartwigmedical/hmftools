package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;

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

    public static synchronized void writeRecord(final BufferedWriter outputWriter, final BlacklistRegion region,
            final AnnotateStatistics stats)
    {
        try
        {
            outputWriter.write(region.getTSVFragment() + '\t' + stats.getTSVFragment());
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

        if(mConfig.BedFile == null)
        {
            MD_LOGGER.error("no BED file specified");
            System.exit(1);
        }

        BufferedWriter outputWriter = createOutputWriter();
        final List<BlacklistRegion> regions =
                BedReader.readFromFile(mConfig.BedFile);

        // TODO(m_cooper): Overlapping
        // TODO(m_cooper): 9
        // TODO(m_cooper): [68999357, 69000000]
        // TODO(m_cooper): [69000000, 69000299]
        // TODO(m_cooper): PosEnd included?

        // TODO(m_cooper): Read entirely contained?

        final ArrayBlockingQueue<BlacklistRegion> jobs = new ArrayBlockingQueue<>(regions.size(), true, regions);
        final List<AnnotateConsumer> annotateConsumers = new ArrayList<>();
        for(int i = 0; i < Math.max(mConfig.Threads, 1); i++)
        {
            final AnnotateConsumer annotateConsumer = new AnnotateConsumer(mConfig, jobs, outputWriter);
            annotateConsumers.add(annotateConsumer);
        }

        TaskExecutor.executeTasks(annotateConsumers, mConfig.Threads);

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
            writer.write(BlacklistRegion.TSV_HEADER + '\t' + AnnotateStatistics.TSV_HEADER);
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
