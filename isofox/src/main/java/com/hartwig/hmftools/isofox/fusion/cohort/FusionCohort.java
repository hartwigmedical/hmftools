package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.ceil;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.PASSING_FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.isofox.cohort.CohortAnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;

public class FusionCohort
{
    private final CohortConfig mConfig;

    private final FusionFilters mFilters;
    private String mFilteredFusionHeader;
    private final ExternalFusionCompare mExternalFusionCompare;

    private final FusionCollection mFusionCollection;

    private BufferedWriter mWriter;

    public static final String PASS_FUSION_FILE_ID = "pass_fusions.csv";

    /* Routines:
        1. generate a cohort file from multiple sample fusion files from Isofox
        2. write a set of filtered/passing fusions for each sample fusion file loaded
        3. run comparison of filtered/passing fusions between Isofox and external fusion files
    */

    public FusionCohort(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mFilters = new FusionFilters(mConfig.Fusions, cmd);
        mFilteredFusionHeader = null;
        mExternalFusionCompare = !mConfig.Fusions.ComparisonSources.isEmpty() ? new ExternalFusionCompare(mConfig) : null;
        mFusionCollection = new FusionCollection(mConfig);
        mWriter = null;

        if(mConfig.Fusions.WriteCombinedFusions)
            intialiseCombinedWriter();
    }

    public void processFusionFiles()
    {
        if(!mConfig.Fusions.GenerateCohort && mConfig.Fusions.ComparisonSources.isEmpty() && !mConfig.Fusions.WriteFilteredFusions)
        {
            ISF_LOGGER.warn("no fusion functions configured");
            return;
        }

        final List<Path> filenames = Lists.newArrayList();

        final CohortAnalysisType fileType = !mConfig.Fusions.ComparisonSources.isEmpty() ? PASSING_FUSION : FUSION;

        if(!formSampleFilenames(mConfig, fileType, filenames))
            return;

        int totalSampleCount = mConfig.SampleData.SampleIds.size();

        int taskId = 0;
        int sampleCount = 0;
        int pairsPerThread = mConfig.Threads > 1 ? (int)ceil(totalSampleCount / (double)mConfig.Threads) : totalSampleCount;

        List<FusionCohortTask> fusionTasks = Lists.newArrayList();

        Map<String,Path> sampleFileMap = Maps.newHashMap();

        for(int i = 0; i < totalSampleCount; ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path fusionFile = filenames.get(i);

            sampleFileMap.put(sampleId, fusionFile);

            ++sampleCount;
            if(sampleCount >= pairsPerThread || i == totalSampleCount - 1)
            {
                fusionTasks.add(new FusionCohortTask(
                        taskId++, mConfig, sampleFileMap,
                        mFilters, mFusionCollection, mExternalFusionCompare, mWriter));

                sampleFileMap = Maps.newHashMap();
                sampleCount = 0;
            }
        }

        ISF_LOGGER.info("loading {} sample fusion files, allocating to {} task(s)", totalSampleCount, fusionTasks.size());

        executeTasks(fusionTasks);

        if(mConfig.Fusions.GenerateCohort)
        {
            mFusionCollection.writeCohortFusions();
            return;
        }

        if(mExternalFusionCompare != null)
            mExternalFusionCompare.close();

        closeBufferedWriter(mWriter);
    }

    private boolean executeTasks(final List<FusionCohortTask> fusionTasks)
    {
        if(mConfig.Threads <= 1)
        {
            fusionTasks.forEach(x -> x.call());
            return true;
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("Isofox-%d").build();

        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
        List<FutureTask> threadTaskList = new ArrayList<FutureTask>();

        for(FusionCohortTask fusionTask : fusionTasks)
        {
            FutureTask futureTask = new FutureTask(fusionTask);

            threadTaskList.add(futureTask);
            executorService.execute(futureTask);
        }

        if(!checkThreadCompletion(threadTaskList))
            return false;

        executorService.shutdown();
        return true;
    }

    private boolean checkThreadCompletion(final List<FutureTask> taskList)
    {
        try
        {
            for (FutureTask futureTask : taskList)
            {
                futureTask.get();
            }
        }
        catch (Exception e)
        {
            ISF_LOGGER.error("task execution error: {}", e.toString());
            e.printStackTrace();
            return false;
        }

        return true;
    }

    private void intialiseCombinedWriter()
    {
        final String outputFile = mConfig.formCohortFilename("combined_fusions.csv");

        try
        {
            mWriter = createBufferedWriter(outputFile, false);
            mWriter.write(String.format("SampleId,%s", mFilteredFusionHeader));
            mWriter.write(",Filter,CohortCount,KnownFusionType");
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
                writer.write(String.format("%s,%s", sampleId, fusion.rawData()));
                writer.write(String.format(",%s,%d,%s",
                        fusion.filter(), fusion.cohortFrequency(), fusion.getKnownFusionType()));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write combined fusion file: {}", e.toString());
        }
    }

}
