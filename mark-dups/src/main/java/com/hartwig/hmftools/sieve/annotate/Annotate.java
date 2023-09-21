package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.NotNull;

public class Annotate
{
    private final AnnotateConfig mConfig;

    public Annotate(final ConfigBuilder configBuilder)
    {
        mConfig = new AnnotateConfig(configBuilder);
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

        final List<AnnotatedBedRecord> bedRecords = BedReader.readFromFile(mConfig.BedFile);

        // TODO(m_cooper): Overlapping
        // TODO(m_cooper): 9
        // TODO(m_cooper): [68999357, 69000000]
        // TODO(m_cooper): [69000000, 69000299]
        // TODO(m_cooper): PosEnd included?

        // TODO(m_cooper): Read entirely contained?

        final List<AnnotateTask> annotateTasks = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            final String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());
            final List<AnnotatedBedRecord> filteredBedRecords = bedRecords.stream().filter(q -> q.getChromosome().equals(stripChrPrefix(chrStr))).collect(Collectors.toList());
            final AnnotateTask annotateTask = new AnnotateTask(chrStr, mConfig, filteredBedRecords);
            annotateTasks.add(annotateTask);
        }

        final List<Callable> callableList = annotateTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        // TODO(m_cooper): Does this wait for all the tasks to finish?
        writeAnnotatedBedFile(bedRecords);

        MD_LOGGER.info("annotate complete");
    }

    private void writeAnnotatedBedFile(@NotNull final List<AnnotatedBedRecord> bedRecords)
    {
        MD_LOGGER.info("Writing output to {}.", mConfig.OutputFile);
        try
        {
            BufferedWriter writer = new BufferedWriter(new FileWriter(mConfig.OutputFile));
            writer.write(AnnotatedBedRecord.TSV_HEADER);
            writer.newLine();
            for(var record : bedRecords)
            {
                writer.write(record.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(Exception e)
        {
            MD_LOGGER.error("An exception was raised while writing the output to {}: {}", mConfig.OutputFile, e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
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
}
