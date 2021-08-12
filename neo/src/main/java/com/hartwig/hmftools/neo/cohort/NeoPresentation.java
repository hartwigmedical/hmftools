package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.EXPRESSION_SCOPE_TRANS;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.cohort.NeoCohortConfig.SAMPLE_TRANS_EXP_FILE;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.bind.BinderConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class NeoPresentation
{
    private final NeoCohortConfig mConfig;
    private final CohortWriters mWriters;
    private final BindScorer mPeptideScorer;
    private final RnaExpressionMatrix mTranscriptExpression;

    public NeoPresentation(final CommandLine cmd)
    {
        mConfig = new NeoCohortConfig(cmd);

        BinderConfig binderConfig = new BinderConfig(cmd);
        mPeptideScorer = new BindScorer(binderConfig);
        mTranscriptExpression = new RnaExpressionMatrix(cmd.getOptionValue(SAMPLE_TRANS_EXP_FILE), EXPRESSION_SCOPE_TRANS);

        mWriters = new CohortWriters(mConfig);
    }

    public void run()
    {
        if(mConfig.SampleIds.isEmpty())
            return;

        if(!mPeptideScorer.loadData())
        {
            NE_LOGGER.error("failed to load scoring data");
            System.exit(1);
        }

        NE_LOGGER.info("processing {} samples", mConfig.SampleIds.size());

        List<NeoSampleBindTask> sampleTasks = Lists.newArrayList();

        for(String sampleId : mConfig.SampleIds)
        {
            NeoSampleBindTask sampleTask = new NeoSampleBindTask(sampleId, mConfig, mPeptideScorer, mTranscriptExpression, mWriters);

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

        mWriters.close();

        NE_LOGGER.info("cohort peptide predictions complete");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        NeoCohortConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        NeoPresentation neoPresentation = new NeoPresentation(cmd);
        neoPresentation.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
