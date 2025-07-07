package com.hartwig.hmftools.isofox.refdata;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.APP_NAME;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class GenerateReferenceData
{
    private final RefDataConfig mConfig;
    private final RefDataWriter mWriter;
    private final EnsemblDataCache mEnsemblDataCache;

    public GenerateReferenceData(final ConfigBuilder configBuilder)
    {
        mConfig = new RefDataConfig(configBuilder);
        mWriter = new RefDataWriter(mConfig);
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
    }

    public void run()
    {
        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            mEnsemblDataCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        }

        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);

        long startTimeMs = System.currentTimeMillis();

        Map<String, List<GeneData>> chrGeneMap = mEnsemblDataCache.getChrGeneDataMap();

        // first execute non-core tasks
        if(mConfig.GenerateGcRatios)
        {
            generateGcRatios(chrGeneMap);
        }

        if(mConfig.GenerateExpectedCounts)
        {
            generateExpectedCounts(chrGeneMap);
        }

        mWriter.close();

        ISF_LOGGER.info("Isofox ref data generation complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private boolean generateExpectedCounts(final Map<String, List<GeneData>> chrGeneMap)
    {
        ISF_LOGGER.info("generating expected transcript counts cache");

        final List<ChrExpectedCountsTask> taskList = Lists.newArrayList();
        final List<Callable<Void>> callableList = Lists.newArrayList();

        for(Map.Entry<String, List<GeneData>> entry : chrGeneMap.entrySet())
        {
            ChrExpectedCountsTask expressionTask = new ChrExpectedCountsTask(mConfig, mEnsemblDataCache, mWriter);
            expressionTask.initialise(entry.getKey(), entry.getValue());
            taskList.add(expressionTask);
            callableList.add(expressionTask);
        }

        return TaskExecutor.executeTasks(callableList, mConfig.Threads);
    }

    private boolean generateGcRatios(final Map<String, List<GeneData>> chrGeneMap)
    {
        ISF_LOGGER.info("generating GC counts cache");

        final List<ExpectedGcRatiosGenerator> taskList = Lists.newArrayList();
        final List<Callable<Void>> callableList = Lists.newArrayList();

        for(Map.Entry<String, List<GeneData>> entry : chrGeneMap.entrySet())
        {
            ExpectedGcRatiosGenerator gcCalcs = new ExpectedGcRatiosGenerator(
                    mConfig, mEnsemblDataCache, entry.getKey(), entry.getValue(), mWriter.getReadGcRatioWriter());

            taskList.add(gcCalcs);
            callableList.add(gcCalcs);
        }

        return TaskExecutor.executeTasks(callableList, mConfig.Threads);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        RefDataConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateReferenceData generateReferenceData = new GenerateReferenceData(configBuilder);
        generateReferenceData.run();
    }
}
