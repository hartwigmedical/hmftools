package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.ceil;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.PASSING_FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.cohort.AnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.isofox.fusion.FusionData;
import com.hartwig.hmftools.isofox.fusion.PassingFusions;

public class FusionCohort
{
    private final CohortConfig mConfig;

    private final PassingFusions mFilters;
    private final BufferedWriter mExternalCompareWriter;
    private final BufferedWriter mUnknownSpliceWriter;

    private final FusionCollection mFusionCollection;

    private BufferedWriter mWriter;

    /* Routines:
        1. generate a cohort file from multiple sample fusion files from Isofox
        2. write a set of filtered/passing fusions for each sample fusion file loaded
        3. run comparison of filtered/passing fusions between Isofox and external fusion files
    */

    public FusionCohort(final CohortConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.loadFromFile(configBuilder);

        mFilters = new PassingFusions(knownFusionCache, mConfig.Fusions.CohortFile);
        mExternalCompareWriter = mConfig.Fusions.ComparisonSource != null ? ExternalFusionCompare.initialiseWriter(mConfig) : null;
        mUnknownSpliceWriter = mConfig.Fusions.FindUnknownSplice ? UnknownSpliceAnalyser.initialiseWriter(mConfig) : null;
        mFusionCollection = new FusionCollection(mConfig);
        mWriter = null;
    }

    public void processFusionFiles()
    {
        if(!mConfig.Fusions.GenerateCohort
        && mConfig.Fusions.ComparisonSource == null
        && !mConfig.Fusions.WriteFilteredFusions
        && !mConfig.Fusions.FindUnknownSplice)
        {
            ISF_LOGGER.warn("no fusion functions configured");
            return;
        }

        final List<Path> filenames = Lists.newArrayList();

        final AnalysisType fileType = mConfig.Fusions.ComparisonSource != null ? PASSING_FUSION : FUSION;

        if(!formSampleFilenames(mConfig, fileType, filenames))
            return;

        int totalSampleCount = mConfig.SampleData.SampleIds.size();

        int taskId = 0;
        int sampleCount = 0;
        int pairsPerThread = mConfig.Threads > 1 ? (int)ceil(totalSampleCount / (double)mConfig.Threads) : totalSampleCount;

        List<FusionCohortTask> fusionTasks = Lists.newArrayList();

        Map<String,Path> sampleFileMap = null;

        for(int i = 0; i < totalSampleCount; ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path fusionFile = filenames.get(i);

            if(mConfig.Fusions.WriteCombinedFusions && mWriter == null)
                intialiseCombinedWriter();

            if(sampleFileMap == null)
                sampleFileMap = Maps.newHashMap();

            sampleFileMap.put(sampleId, fusionFile);

            ++sampleCount;
            if(sampleCount >= pairsPerThread || i == totalSampleCount - 1)
            {
                fusionTasks.add(new FusionCohortTask(
                        taskId++, mConfig, sampleFileMap, mFilters, mFusionCollection,
                        mWriter, mExternalCompareWriter, mUnknownSpliceWriter));

                sampleFileMap = null;
                sampleCount = 0;
            }
        }

        ISF_LOGGER.info("loading ({}) sample fusion files, allocating to {} task(s)", totalSampleCount, fusionTasks.size());

        final List<Callable> callableList = fusionTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        if(mConfig.Fusions.GenerateCohort)
        {
            mFusionCollection.writeCohortFusions();
            return;
        }

        closeBufferedWriter(mWriter);

        if(mExternalCompareWriter != null)
        {
            fusionTasks.forEach(x -> ExternalFusionCompare.writeResults(mExternalCompareWriter, x.getExternalCompare().getResults()));
            closeBufferedWriter(mExternalCompareWriter);
        }

        closeBufferedWriter(mUnknownSpliceWriter);

        ISF_LOGGER.info("fusion cohort analysis complete");
    }

    private void intialiseCombinedWriter()
    {
        final String outputFile = mConfig.formCohortFilename("combined_fusions.csv");

        try
        {
            mWriter = createBufferedWriter(outputFile, false);
            mWriter.write(String.format("SampleId,%s", FusionData.csvHeader(true)));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialise combined fusion file({}): {}", outputFile, e.toString());
        }
    }

    public synchronized static void writeCombinedFusions(final BufferedWriter writer, final String sampleId, final List<FusionData> fusions)
    {
        if(writer == null)
            return;

        try
        {
            for (FusionData fusion : fusions)
            {
                writer.write(String.format("%s,%s", sampleId, fusion.toCsv(true)));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write combined fusion file: {}", e.toString());
        }
    }

}
