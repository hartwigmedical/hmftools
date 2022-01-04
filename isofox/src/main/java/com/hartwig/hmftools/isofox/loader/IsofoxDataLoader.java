package com.hartwig.hmftools.isofox.loader;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.loader.DataLoadType.GENE_EXPRESSION;
import static com.hartwig.hmftools.isofox.loader.DataLoadType.NOVEL_JUNCTION;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.isofox.expression.cohort.CohortGenePercentiles;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class IsofoxDataLoader
{
    private final DataLoaderConfig mConfig;
    private final AltSjCohortCache mAltSjCohortCache;
    private final CohortGenePercentiles mGeneDistribution;

    public IsofoxDataLoader(final CommandLine cmd)
    {
        mConfig = new DataLoaderConfig(cmd);

        if(mConfig.loadDataType(NOVEL_JUNCTION) && mConfig.AltSjCohortFile != null)
            mAltSjCohortCache = new AltSjCohortCache(mConfig.AltSjCohortFile);
        else
            mAltSjCohortCache = null;

        if(mConfig.loadDataType(GENE_EXPRESSION) && mConfig.GeneDistributionFile != null)
            mGeneDistribution = new CohortGenePercentiles(mConfig.GeneDistributionFile);
        else
            mGeneDistribution = null;
    }

    public boolean load()
    {
        if(mConfig.SampleIds.isEmpty())
        {
            ISF_LOGGER.error("no sample IDs configured");
            return false;
        }

        if(mConfig.DbAccess == null)
        {
            ISF_LOGGER.error("invalid DB connection");
            return false;
        }

        List<SampleLoaderTask> sampleTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            for(int i = 0; i < min(mConfig.SampleIds.size(), mConfig.Threads); ++i)
            {
                sampleTasks.add(new SampleLoaderTask(i, mConfig, mAltSjCohortCache, mGeneDistribution));
            }

            int taskIndex = 0;
            for(String sampleId : mConfig.SampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            SampleLoaderTask sampleTask = new SampleLoaderTask(0, mConfig, mAltSjCohortCache, mGeneDistribution);

            sampleTask.getSampleIds().addAll(mConfig.SampleIds);

            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = DataLoaderConfig.createCmdLineOptions();
        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        setLogLevel(cmd);

        IsofoxDataLoader dataLoader = new IsofoxDataLoader(cmd);

        if(!dataLoader.load())
        {
            ISF_LOGGER.info("Isofox data loading failed");
            return;
        }

        ISF_LOGGER.info("Isofox data loading complete");
    }
}
