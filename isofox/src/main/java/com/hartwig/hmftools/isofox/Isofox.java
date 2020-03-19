package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;

import static com.hartwig.hmftools.isofox.ChromosomeGeneTask.PERF_FIT;
import static com.hartwig.hmftools.isofox.ChromosomeGeneTask.PERF_FRAG_LENGTH;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_LEVEL;
import static com.hartwig.hmftools.isofox.IsofoxConfig.createCmdLineOptions;
import static com.hartwig.hmftools.isofox.TaskType.FRAGMENT_LENGTHS;
import static com.hartwig.hmftools.isofox.TaskType.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.common.FragmentSizeCalcs.setConfigFragmentLengthData;
import static com.hartwig.hmftools.isofox.results.SummaryStats.createSummaryStats;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengths;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.isofox.common.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedCountsCache;
import com.hartwig.hmftools.isofox.gc.GcBiasAdjuster;
import com.hartwig.hmftools.isofox.gc.GcRatioCounts;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.isofox.results.SummaryStats;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class Isofox
{
    private final IsofoxConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final EnsemblDataCache mGeneTransCache;
    private final GcBiasAdjuster mGcBiasAdjuster;
    private final ExpectedCountsCache mExpectedCountsCache;
    private final ExecutorService mExecutorService;

    private final List<int[]> mFragmentLengthDistribution;

    public Isofox(final IsofoxConfig config, final CommandLine cmd)
    {
        mConfig = config;

        mGcBiasAdjuster = new GcBiasAdjuster(mConfig);

        mResultsWriter = new ResultsWriter(mConfig);

        mGeneTransCache = new EnsemblDataCache(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), RefGenomeVersion.HG37);

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        }

        mGeneTransCache.setRequiredData(true, false, false, mConfig.CanonicalTranscriptOnly);
        mGeneTransCache.load(false);

        if(mConfig.Threads > 1)
        {
            final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("RnaExp-%d").build();
            mExecutorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
        }
        else
        {
            mExecutorService = null;
        }

        mExpectedCountsCache = mConfig.ExpCountsFile != null ? new ExpectedCountsCache(mConfig) : null;

        mFragmentLengthDistribution = Lists.newArrayList();
    }

    public void runAnalysis()
    {
        ISF_LOGGER.info("sample({}) running RNA expression analysis", mConfig.SampleId);

        if(mGcBiasAdjuster.enabled())
        {
            mGcBiasAdjuster.loadData();
            mGcBiasAdjuster.generateDepthCounts(mGeneTransCache.getChrGeneDataMap());
            return; // for now
        }

        // allocate work at the chromosome level
        List<ChromosomeGeneTask> chrTasks = Lists.newArrayList();

        for(Map.Entry<String,List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final List<EnsemblGeneData> geneDataList = entry.getValue();

            final String chromosome = entry.getKey();

            if (mConfig.skipChromosome(chromosome) || geneDataList.isEmpty())
                continue;

            ChromosomeGeneTask chrGeneTask = new ChromosomeGeneTask(
                    mConfig, chromosome, geneDataList, mGeneTransCache, mResultsWriter, mExpectedCountsCache);

            chrTasks.add(chrGeneTask);
        }

        if(mConfig.requireFragmentLengthCalcs())
        {
            calcFragmentLengths(chrTasks);

            if(mConfig.WriteFragmentLengthsOnly)
            {
                mResultsWriter.close();
                return;
            }
        }

        boolean validExecution = executeChromosomeTask(chrTasks, TRANSCRIPT_COUNTS);

        if(!validExecution)
            return;

        // final reporting
        if(!mConfig.generateExpRatesOnly())
        {
            int totalReadsProcessed = chrTasks.stream().mapToInt(x -> x.getBamReader().totalBamCount()).sum();
            ISF_LOGGER.info("read {} total BAM records", totalReadsProcessed);

            int totalFragCount = chrTasks.stream().mapToInt(x -> x.getTotalFragmentCount()).sum();
            int enrichedGeneFragCount = chrTasks.stream().mapToInt(x -> x.getEnrichedGenesFragmentCount()).sum();

            GcRatioCounts nonEnrichedGcRatioCounts = new GcRatioCounts();
            chrTasks.forEach(x -> nonEnrichedGcRatioCounts.mergeRatioCounts(x.getNonEnrichedGcRatioCounts().getRatioCounts()));
            double medianGCRatio = nonEnrichedGcRatioCounts.getPercentileRatio(0.5);

            final SummaryStats summaryStats = createSummaryStats(
                    totalFragCount, enrichedGeneFragCount, medianGCRatio, mFragmentLengthDistribution, mConfig.ReadLength);

            mResultsWriter.writeSummaryStats(summaryStats);
        }

        if(mConfig.WriteReadGcRatios)
        {
            GcRatioCounts ratioCounts = new GcRatioCounts();
            chrTasks.forEach(x -> ratioCounts.mergeRatioCounts(x.getBamReader().getGcRatioCounts().getRatioCounts()));

            GcRatioCounts.writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), null, ratioCounts.getRatioCounts());
        }

        mResultsWriter.close();

        final PerformanceCounter[] perfCounters = chrTasks.get(0).getPerfCounters();

        for(int i = 1; i < chrTasks.size(); ++i)
        {
            final PerformanceCounter[] chrPCs = chrTasks.get(i).getPerfCounters();

            for(int j = 0; j < perfCounters.length; ++j)
            {
                perfCounters[j].merge(chrPCs[j]);
            }
        }

        Arrays.stream(perfCounters).forEach(x -> x.logStats());

        if(mConfig.RunPerfChecks)
        {
            // log 10 slowest times and their interval names
            final List<Double> fitTimes = perfCounters[PERF_FIT].getTimes();
            final List<String> fitGenes = perfCounters[PERF_FIT].getTimeNames();

            if(fitTimes.size() >= 10 && fitGenes.size() == fitTimes.size())

            for (int i = fitTimes.size() - 1; i >= fitTimes.size() - 10; --i)
            {
                ISF_LOGGER.info(String.format("fit times: geneSet(%s) time(%.3f)", fitGenes.get(i), fitTimes.get(i)));
            }
        }
    }

    private boolean executeChromosomeTask(final List<ChromosomeGeneTask> chrTasks, TaskType taskType)
    {
        chrTasks.forEach(x -> x.setTaskType(taskType));

        if(mConfig.Threads <= 1 || mExecutorService == null)
        {
            chrTasks.forEach(x -> x.call());
            return true;
        }

        List<FutureTask> threadTaskList = new ArrayList<FutureTask>();

        for(ChromosomeGeneTask chrGeneTask : chrTasks)
        {
                FutureTask futureTask = new FutureTask(chrGeneTask);

                threadTaskList.add(futureTask);
                mExecutorService.execute(futureTask);
        }

        return checkThreadCompletion(threadTaskList);
    }

    private void calcFragmentLengths(final List<ChromosomeGeneTask> chrTasks)
    {
        // for now a way of only calculating fragment lengths and nothing more
        boolean validExecution = executeChromosomeTask(chrTasks, FRAGMENT_LENGTHS);

        if(!validExecution)
            return;

        // merge results from all chromosomes
        int maxReadLength = 0;
        for(final ChromosomeGeneTask chrGeneTask : chrTasks)
        {
            final FragmentSizeCalcs fragSizeCalcs = chrGeneTask.getFragSizeCalcs();
            maxReadLength = max(maxReadLength, fragSizeCalcs.getMaxReadLength());
            FragmentSizeCalcs.mergeData(mFragmentLengthDistribution, fragSizeCalcs);
        }

        if(mConfig.UseCalculatedFragmentLengths)
            setConfigFragmentLengthData(mConfig, maxReadLength, mFragmentLengthDistribution);

        if (mConfig.WriteFragmentLengths && !mConfig.FragmentLengthsByGene)
        {
            FragmentSizeCalcs.writeFragmentLengths(mResultsWriter.getFragmentLengthWriter(), mFragmentLengthDistribution, null);
        }

        if(mConfig.WriteFragmentLengthsOnly)
        {
            final PerformanceCounter perfCounter = chrTasks.get(0).getPerfCounters()[PERF_FRAG_LENGTH];

            for(int i = 1; i < chrTasks.size(); ++i)
            {
                perfCounter.merge(chrTasks.get(i).getPerfCounters()[PERF_FRAG_LENGTH]);
            }

            perfCounter.logStats();
        }
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

        mExecutorService.shutdown();
        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final VersionInfo version = new VersionInfo("isofox.version");
        ISF_LOGGER.info("Isofox version: {}", version.version());

        final Options options = createCmdLineOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }
        else if(cmd.hasOption(LOG_LEVEL))
        {
            Configurator.setRootLevel(Level.valueOf(cmd.getOptionValue(LOG_LEVEL)));
        }

        IsofoxConfig config = new IsofoxConfig(cmd);

        if(!config.isValid())
        {
            ISF_LOGGER.error("missing config options, exiting");
            return;
        }

        Isofox isofox = new Isofox(config, cmd);
        isofox.runAnalysis();

        ISF_LOGGER.info("RNA expression analysis complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
