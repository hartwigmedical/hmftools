package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.createCmdLineOptions;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class RnaExpression
{
    private final RnaExpConfig mConfig;
    private final RnaBamReader mRnaBamReader;
    private final ResultsWriter mResultsWriter;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final GcBiasAdjuster mGcBiasAdjuster;
    private final FragmentSizeCalcs mFragmentSizeCalcs;
    private final ExecutorService mExecutorService;

    public RnaExpression(final CommandLine cmd)
    {
        mConfig = new RnaExpConfig(cmd);

        mRnaBamReader = new RnaBamReader(mConfig);
        mGcBiasAdjuster = new GcBiasAdjuster(mConfig);

        mResultsWriter = new ResultsWriter(mConfig);

        mGeneTransCache = new SvGeneTranscriptCollection();
        mGeneTransCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        }

        mGeneTransCache.setRequiredData(true, false, false, mConfig.CanonicalTranscriptOnly);
        mGeneTransCache.loadEnsemblData(false);

        mFragmentSizeCalcs = new FragmentSizeCalcs(mConfig, mGeneTransCache);

        if(mConfig.Threads > 1)
        {
            final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("RnaExp-%d").build();
            mExecutorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
        }
        else
        {
            mExecutorService = null;
        }
    }

    public void runAnalysis()
    {
        RE_LOGGER.info("sample({}) running RNA expression analysis", mConfig.SampleId);

        if(mGcBiasAdjuster.enabled())
        {
            mGcBiasAdjuster.loadData();
            mGcBiasAdjuster.generateDepthCounts(mRnaBamReader, mGeneTransCache.getChrGeneDataMap());
            return; // for now
        }

        if(mConfig.FragmentLengthMinCount > 0)
        {
            mFragmentSizeCalcs.calcSampleFragmentSize();

            if (mConfig.WriteFragmentLengths)
                return;
        }

        // prepared writers for threaded execution
        BufferedWriter expRatesWriter = mResultsWriter.initialiseExpRatesWriter();
        BufferedWriter readsWriter = null; // mResultsWriter.initialiseExpRatesWriter();

        List<FutureTask> taskList = new ArrayList<FutureTask>();

        for(Map.Entry<String,List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final List<EnsemblGeneData> geneDataList = entry.getValue();

            final String chromosome = entry.getKey();

            if(mConfig.skipChromosome(chromosome))
                continue;

            ChromosomeGeneTask chrGeneTask = new ChromosomeGeneTask(mConfig, chromosome, geneDataList, mGeneTransCache, mResultsWriter);
            chrGeneTask.setWriters(expRatesWriter, readsWriter);

            if(mExecutorService != null)
            {
                FutureTask futureTask = new FutureTask(chrGeneTask);

                taskList.add(futureTask);
                mExecutorService.execute(futureTask);
            }
            else
            {
                chrGeneTask.analyseGenes();
            }
        }

        if(mExecutorService != null)
        {
            // Wait until all results are available and combine them at the same time
            for (FutureTask futureTask : taskList)
            {
                try
                {
                    futureTask.get();
                }
                catch (Exception e)
                {
                    RE_LOGGER.error("task execution error: {}", e.toString());
                }
            }

            mExecutorService.shutdown();
        }

        mResultsWriter.close();
        mRnaBamReader.close();
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createCmdLineOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        if(!RnaExpConfig.checkValid(cmd))
        {
            RE_LOGGER.error("missing config options, exiting");
            return;
        }

        RnaExpression rnaExpression = new RnaExpression(cmd);
        rnaExpression.runAnalysis();

        RE_LOGGER.info("RNA expression analysis complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
