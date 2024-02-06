package com.hartwig.hmftools.neo.score;

import static java.lang.Math.min;

import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class NeoScorer
{
    private final NeoScorerConfig mConfig;
    private final NeoDataWriter mWriters;
    private final ReferenceData mReferenceData;

    public NeoScorer(final ConfigBuilder configBuilder)
    {
        mConfig = new NeoScorerConfig(configBuilder);
        mReferenceData = new ReferenceData(configBuilder);

        mWriters = new NeoDataWriter(mConfig);
    }

    public void run()
    {
        if(mConfig.Samples.isEmpty())
            return;

        if(!mReferenceData.PeptideScorer.loadScoringData())
        {
            NE_LOGGER.error("failed to load scoring data");
            System.exit(1);
        }

        NE_LOGGER.info("running neoepitope scoring for {}",
                mConfig.Samples.size() == 1 ? mConfig.Samples.get(0).TumorId : String.format("%d samples", mConfig.Samples.size()));

        List<NeoScorerTask> sampleTasks = Lists.newArrayList();

        for(int i = 0; i < min(mConfig.Threads, mConfig.Samples.size()); ++i)
        {
            sampleTasks.add(new NeoScorerTask(i, mConfig, mReferenceData, mWriters));
        }

        int taskIndex = 0;
        for(SampleData sampleData : mConfig.Samples)
        {
            sampleTasks.get(taskIndex).addSample(sampleData);

            ++taskIndex;
            if(taskIndex == sampleTasks.size())
                taskIndex = 0;
        }

        final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        mWriters.close();

        NE_LOGGER.info("Neo peptide scoring complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        NeoScorerConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        NeoScorer neoScorer = new NeoScorer(configBuilder);
        neoScorer.run();
    }
}
