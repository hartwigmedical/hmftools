package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class NeoEpitopeFinder
{
    private final NeoConfig mConfig;

    private final EnsemblDataCache mGeneTransCache;

    public NeoEpitopeFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new NeoConfig(configBuilder);

        mGeneTransCache = new EnsemblDataCache(configBuilder);
        mGeneTransCache.setRequiredData(true, false, false, false);
        mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        mGeneTransCache.load(false);
        mGeneTransCache.createGeneNameIdMap();
    }

    public void run()
    {
        if(mConfig.Samples.isEmpty())
            return;

        if(mConfig.Samples.size() == 1)
        {
            NE_LOGGER.info("processing sample({})", mConfig.Samples.get(0));
        }
        else
        {
            NE_LOGGER.info("processing {} samples", mConfig.Samples.size());
        }

        // check required inputs and config
        List<NeoSampleTask> sampleTasks = Lists.newArrayList();

        for(final String sample : mConfig.Samples)
        {
            NeoSampleTask sampleTask = new NeoSampleTask(sample, mConfig, mGeneTransCache);

            sampleTasks.add(sampleTask);
        }

        if(mConfig.Threads > 1)
        {
            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            sampleTasks.forEach(x -> x.processSample());
        }

        NE_LOGGER.info("Neo-epitope annotations complete");
    }

    public static BufferedWriter initialiseNeoepitopeWriter(final String outputDir, final String sampleId)
    {
        if(outputDir.isEmpty())
            return null;

        try
        {
            String outputFileName = sampleId != null ?
                    NeoEpitopeFile.generateFilename(outputDir, sampleId) : outputDir + "IMU_NEO_EPITOPES.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(sampleId == null)
                writer.write("SampleId,");

            writer.write(NeoEpitopeFile.header());
            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            NE_LOGGER.error("error initialising neo-epitope output file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeNeoepitopes(
            final BufferedWriter writer, final String sampleId, boolean isCohort,
            int neId, final NeoEpitope neData, final Set<String> upTransNames, final Set<String> downTransNames)
    {
        if(writer == null)
            return;

        try
        {
            if(isCohort)
                writer.write(String.format("%s,", sampleId));

            final NeoEpitopeFile neFile = neData.toFile(neId, upTransNames, downTransNames);
            writer.write(NeoEpitopeFile.toString(neFile));
            writer.newLine();
        }
        catch (final IOException e)
        {
            NE_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        NeoConfig.addCmdLineArgs(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        NeoEpitopeFinder neoEpitopeFinder = new NeoEpitopeFinder(configBuilder);
        neoEpitopeFinder.run();
    }
}
