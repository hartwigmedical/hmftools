package com.hartwig.hmftools.neo.scorer;

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

//
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

        NE_LOGGER.info("loading cohort transcript expression");

        RnaExpressionMatrix transcriptExpression = mConfig.CohortSampleTpmFile != null ?
                new RnaExpressionMatrix(mConfig.CohortSampleTpmFile, EXPRESSION_SCOPE_TRANS) : null;

        TpmMediansCache tpmMediansCache = new TpmMediansCache(mConfig.CohortTpmMediansFile);

        NE_LOGGER.info("running neoepitope scoring for {}",
                mConfig.Samples.size() == 1 ? mConfig.Samples.get(0).Id : String.format("%d samples", mConfig.Samples.size()));

        List<NeoScorerTask> sampleTasks = Lists.newArrayList();

        for(SampleData sampleData : mConfig.Samples)
        {
            NeoScorerTask sampleTask = new NeoScorerTask(sampleData, mConfig, mPeptideScorer, transcriptExpression, tpmMediansCache, mWriters);

            sampleTasks.add(sampleTask);
        }

        if(mConfig.Threads > 1)
        {
            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            sampleTasks.forEach(x -> x.call());
        }

        mWriters.close();

        NE_LOGGER.info("cohort peptide predictions complete");
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
