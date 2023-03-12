package com.hartwig.hmftools.neo.scorer;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.EXPRESSION_SCOPE_TRANS;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.bind.ScoreConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class NeoScorer
{
    private final NeoScorerConfig mConfig;
    private final NeoDataWriter mWriters;
    private final BindScorer mPeptideScorer;

    public NeoScorer(final CommandLine cmd)
    {
        mConfig = new NeoScorerConfig(cmd);

        mPeptideScorer = new BindScorer(new ScoreConfig(cmd));

        mWriters = new NeoDataWriter(mConfig);
    }

    public void run()
    {
        if(mConfig.Samples.isEmpty())
            return;

        if(!mPeptideScorer.loadScoringData())
        {
            NE_LOGGER.error("failed to load scoring data");
            System.exit(1);
        }

        RnaExpressionMatrix transcriptExpression = null;

        if(mConfig.CohortSampleTpmFile != null)
        {
            NE_LOGGER.info("loading cohort transcript expression");
            transcriptExpression = new RnaExpressionMatrix(mConfig.CohortSampleTpmFile, EXPRESSION_SCOPE_TRANS);
        }

        NE_LOGGER.info("loading cohort transcript medians");

        TpmMediansCache tpmMediansCache = new TpmMediansCache(mConfig.CohortTpmMediansFile);

        NE_LOGGER.info("running neoepitope scoring for {}",
                mConfig.Samples.size() == 1 ? mConfig.Samples.get(0).Id : String.format("%d samples", mConfig.Samples.size()));

        List<NeoScorerTask> sampleTasks = Lists.newArrayList();

        for(int i = 0; i < min(mConfig.Threads, mConfig.Samples.size()); ++i)
        {
            sampleTasks.add(new NeoScorerTask(mConfig, mPeptideScorer, transcriptExpression, tpmMediansCache, mWriters));
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
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        mWriters.close();

        NE_LOGGER.info("cohort neoepitope peptide scoring complete");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        NeoScorerConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        NeoScorer neoScorer = new NeoScorer(cmd);
        neoScorer.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
