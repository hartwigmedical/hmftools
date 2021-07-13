package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sigs.SigUtils.convertToPercentages;
import static com.hartwig.hmftools.common.sigs.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.BamFragmentReader.PERF_FIT;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.createCmdLineOptions;
import static com.hartwig.hmftools.isofox.IsofoxConfig.validConfigPaths;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_GC_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_TRANS_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.NEO_EPITOPES;
import static com.hartwig.hmftools.isofox.IsofoxFunction.READ_COUNTS;
import static com.hartwig.hmftools.isofox.TaskType.APPLY_GC_ADJUSTMENT;
import static com.hartwig.hmftools.isofox.TaskType.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs.setConfigFragmentLengthData;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.writeReadGcRatioCounts;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.expression.TranscriptExpression.calcTpmFactors;
import static com.hartwig.hmftools.isofox.expression.TranscriptExpression.setTranscriptsPerMillion;
import static com.hartwig.hmftools.isofox.results.SummaryStats.createSummaryStats;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.adjusts.GcTranscriptCalculator;
import com.hartwig.hmftools.isofox.common.BamReadCounter;
import com.hartwig.hmftools.isofox.common.FragmentType;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.isofox.expression.ExpectedCountsCache;
import com.hartwig.hmftools.isofox.expression.ExpressionCacheTask;
import com.hartwig.hmftools.isofox.expression.GeneCollectionSummary;
import com.hartwig.hmftools.isofox.fusion.ChimericReadCache;
import com.hartwig.hmftools.isofox.fusion.ChimericStats;
import com.hartwig.hmftools.isofox.fusion.FusionTaskManager;
import com.hartwig.hmftools.isofox.neo.NeoEpitopeReader;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class Isofox
{
    private final IsofoxConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final EnsemblDataCache mGeneTransCache;
    private final ExpectedCountsCache mExpectedCountsCache;
    private final GcTranscriptCalculator mGcTranscriptCalcs;
    private final FusionTaskManager mFusionTaskManager;

    private int mMaxObservedReadLength;
    private final List<FragmentSize> mFragmentLengthDistribution;

    public Isofox(final IsofoxConfig config, final CommandLine cmd)
    {
        mConfig = config;

        mResultsWriter = new ResultsWriter(mConfig);

        mGeneTransCache = new EnsemblDataCache(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), config.RefGenVersion);

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        }

        mGeneTransCache.setRequiredData(true, false, false, mConfig.CanonicalTranscriptOnly);
        mGeneTransCache.load(false);

        mExpectedCountsCache = mConfig.ExpCountsFile != null || mConfig.ApplyGcBiasAdjust ? new ExpectedCountsCache(mConfig) : null;

        mGcTranscriptCalcs = mConfig.runFunction(EXPECTED_GC_COUNTS) || mConfig.ApplyGcBiasAdjust ?
                new GcTranscriptCalculator(mConfig, mGeneTransCache) : null;

        mFusionTaskManager = mConfig.runFunction(FUSIONS) ? new FusionTaskManager(mConfig, mGeneTransCache) : null;

        mMaxObservedReadLength = 0;
        mFragmentLengthDistribution = Lists.newArrayList();
    }

    public boolean runAnalysis()
    {
        if(mConfig.runFunction(FUSIONS) && mConfig.Fusions.ChimericReadsFile != null)
        {
            mFusionTaskManager.processCachedFragments(ChimericReadCache.loadChimericReads(mConfig.Fusions.ChimericReadsFile));
            return true;
        }

        // all other routines split work by chromosome
        final Map<String,List<GeneData>> chrGeneMap = getChromosomeGeneLists();

        if(chrGeneMap.isEmpty())
        {
            ISF_LOGGER.error("no chromosome tasks created");
            return false;
        }

        // first execute non-core tasks
        if(mConfig.runFunction(EXPECTED_GC_COUNTS))
        {
            boolean status = generateGcRatios(chrGeneMap);
            return status;
        }

        if(mConfig.runFunction(EXPECTED_TRANS_COUNTS))
        {
            boolean status = generateExpectedCounts(chrGeneMap);
            mResultsWriter.close();
            return status;
        }

        if(mConfig.runFunction(NEO_EPITOPES))
        {
            NeoEpitopeReader neReader = new NeoEpitopeReader(mConfig, mGeneTransCache);
            neReader.calcFragmentSupport();
            return true;
        }

        if(mConfig.runFunction(READ_COUNTS))
        {
            return countBamReads(chrGeneMap);
        }

        // BAM processing for the key routines - novel junctions, fusions and gene expression
        if(!allocateBamFragments(chrGeneMap))
            return false;

        return true;
    }

    private boolean allocateBamFragments(final Map<String,List<GeneData>> chrGeneMap)
    {
        ISF_LOGGER.info("sample({}) running RNA analysis", mConfig.SampleId);

        if(mExpectedCountsCache != null && !mExpectedCountsCache.isValid())
        {
            ISF_LOGGER.warn("invalid expected counts cache");
            return false;
        }

        if(mConfig.requireFragmentLengthCalcs())
        {
            calcFragmentLengths(chrGeneMap);

            if(mConfig.WriteFragmentLengthsByGene)
            {
                mResultsWriter.close();
                return true;
            }
        }

        final List<BamFragmentReader> chrTasks = Lists.newArrayList();
        final List<Callable> callableList = Lists.newArrayList();
        final List<String> chromosomes = Lists.newArrayList();

        // process any enriched genes first, then add the rest in order of decreasing length
        chromosomes.add(mConfig.RefGenVersion.versionedChromosome("14"));
        chromosomes.add(mConfig.RefGenVersion.versionedChromosome("3"));
        chromosomes.add(mConfig.RefGenVersion.versionedChromosome("6"));
        chromosomes.add(mConfig.RefGenVersion.versionedChromosome("9"));

        Arrays.stream(HumanChromosome.values())
                .map(chromosome -> mConfig.RefGenVersion.versionedChromosome(chromosome.toString()))
                .filter(chromosome -> !chromosomes.contains(chromosome))
                .forEach(chromosome -> chromosomes.add(chromosome));

        for(String chromosome : chromosomes)
        {
            List<GeneData> geneDataList = chrGeneMap.get(chromosome);

            if(geneDataList == null)
                continue;

            BamFragmentReader bamReaderTask = new BamFragmentReader(
                    mConfig, chromosome, geneDataList, mGeneTransCache, mResultsWriter,
                    mFusionTaskManager, mExpectedCountsCache, mGcTranscriptCalcs);

            chrTasks.add(bamReaderTask);
            callableList.add(bamReaderTask);
        }

        chrTasks.forEach(x -> x.setTaskType(TRANSCRIPT_COUNTS));

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            return false;

        int totalReadsProcessed = chrTasks.stream().mapToInt(x -> x.totalReadCount()).sum();
        ISF_LOGGER.info("read {} total BAM records", totalReadsProcessed);

        if(!mConfig.runFusionsOnly())
        {
            // post processing for summary stats and gene expression data
            processBamFragments(chrTasks, callableList);
        }

        if(mConfig.runFunction(FUSIONS))
        {
            // extract all chimeric reads and associated data for fusion calling
            ChimericStats chimericStats = new ChimericStats();
            chrTasks.forEach(x -> chimericStats.merge(x.getChimericStats()));
            ISF_LOGGER.info("overall chimeric stats: {}", chimericStats);
            mFusionTaskManager.close();
        }

        final List<PerformanceCounter[]> perfCounters = chrTasks.stream().map(x -> x.getPerfCounters()).collect(Collectors.toList());
        chrTasks.clear();

        logPerformanceStats(perfCounters);
        return true;
    }

    private void processBamFragments(final List<BamFragmentReader> chrTasks, final List<Callable> callableList)
    {
        int[] totalCounts = new int[typeAsInt(FragmentType.MAX)];

        for (int i = 0; i < totalCounts.length; ++i)
        {
            final int fragIndex = i;
            totalCounts[i] += chrTasks.stream().mapToInt(x -> x.getCombinedCounts()[fragIndex]).sum();
        }

        int enrichedGeneFragCount = chrTasks.stream().mapToInt(x -> x.getEnrichedGenesFragmentCount()).sum();

        GcRatioCounts nonEnrichedGcRatioCounts = new GcRatioCounts();
        chrTasks.forEach(x -> nonEnrichedGcRatioCounts.mergeRatioCounts(x.getNonEnrichedGcRatioCounts().getCounts()));
        double medianGCRatio = nonEnrichedGcRatioCounts.getPercentileRatio(0.5);

        if (mConfig.ApplyGcBiasAdjust)
        {
            applyGcAdjustments(chrTasks, callableList, nonEnrichedGcRatioCounts);
        }

        if(mConfig.runFunction(IsofoxFunction.TRANSCRIPT_COUNTS))
        {
            final RnaStatistics summaryStats = createSummaryStats(
                    totalCounts, enrichedGeneFragCount,
                    medianGCRatio, mFragmentLengthDistribution, mMaxObservedReadLength > 0 ? mMaxObservedReadLength : mConfig.ReadLength);

            mResultsWriter.writeSummaryStats(summaryStats);
        }

        if (mConfig.WriteGcData)
        {
            GcRatioCounts combinedGcRatioCounts = new GcRatioCounts();
            chrTasks.forEach(x -> combinedGcRatioCounts.mergeRatioCounts(x.getGcRatioCounts().getCounts()));

            writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "ALL", combinedGcRatioCounts.getCounts(), false);
            double[] percentData = new double[combinedGcRatioCounts.size()];

            copyVector(combinedGcRatioCounts.getCounts(), percentData);
            convertToPercentages(percentData);
            writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "ALL_PERC", percentData, true);

            if (!mConfig.EnrichedGeneIds.isEmpty())
                writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "NON_ENRICHED", nonEnrichedGcRatioCounts.getCounts(), false);

            if (mConfig.ApplyGcBiasAdjust)
            {
                writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "TRANS_FIT_EXPECTED",
                        mGcTranscriptCalcs.getTranscriptFitGcCounts().getCounts(), false);

                copyVector(mGcTranscriptCalcs.getTranscriptFitGcCounts().getCounts(), percentData);
                convertToPercentages(percentData);
                writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "TRANS_FIT_EXPECTED_PERC", percentData, true);

                writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "ADJUSTMENTS",
                        mGcTranscriptCalcs.getGcRatioAdjustments(), true);
            }
        }

        // calculate a TPM for all transcripts before results are written
        final List<GeneCollectionSummary> geneSummaryData = Lists.newArrayList();
        chrTasks.stream().forEach(x -> geneSummaryData.addAll(x.getGeneCollectionSummaryData()));

        double[] tpmFactors = calcTpmFactors(geneSummaryData, mConfig.EnrichedGeneIds);
        chrTasks.forEach(x -> setTranscriptsPerMillion(x.getGeneCollectionSummaryData(), tpmFactors));

        chrTasks.forEach(x -> x.writeResults());
        mResultsWriter.close();
    }

    private void applyGcAdjustments(
            final List<BamFragmentReader> chrTasks, final List<Callable> callableList, final GcRatioCounts actualGcCounts)
    {
        ISF_LOGGER.info("applying GC adjustments and transcript re-fit");

        // not thread safe at the moment
        chrTasks.forEach(x -> x.applyGcAdjustment());

        ISF_LOGGER.debug("total({}) transcript expected GC counts from fit", String.format("%.0f",
                mGcTranscriptCalcs.getTranscriptFitGcCounts().getCountsTotal()));

        // gather up global expected counts
        mGcTranscriptCalcs.calcGcRatioAdjustments(actualGcCounts);

        // now re-fit all transcripts
        chrTasks.forEach(x -> x.setTaskType(APPLY_GC_ADJUSTMENT));
        TaskExecutor.executeTasks(callableList, mConfig.Threads);
    }

    private Map<String,List<GeneData>> getChromosomeGeneLists()
    {
        if(mConfig.SpecificRegions.isEmpty() && mConfig.SpecificChromosomes.isEmpty())
            return mGeneTransCache.getChrGeneDataMap();

        final Map<String,List<GeneData>> chrGeneMap = Maps.newHashMap();

        if(mConfig.SpecificRegions.isEmpty())
        {
            mGeneTransCache.getChrGeneDataMap().entrySet().stream()
                    .filter(x -> mConfig.SpecificChromosomes.contains(x.getKey()))
                    .filter(x -> !x.getValue().isEmpty())
                    .forEach(x -> chrGeneMap.put(x.getKey(), x.getValue()));
        }
        else
        {
            for(final BaseRegion region : mConfig.SpecificRegions)
            {
                List<GeneData> geneDataList = mGeneTransCache.getChrGeneDataMap().get(region.Chromosome);

                List<GeneData> regionGeneList = geneDataList.stream()
                        .filter(x -> positionsOverlap(region.start(), region.end(), x.GeneStart, x.GeneEnd))
                        .collect(Collectors.toList());

                chrGeneMap.put(region.Chromosome, regionGeneList);
            }
        }

        return chrGeneMap;
    }

    private void calcFragmentLengths(final Map<String,List<GeneData>> chrGeneMap)
    {
        int requiredFragCount = mConfig.FragmentLengthSamplingCount / chrGeneMap.size(); // split evenly amongst chromosomes

        final List<FragmentSizeCalcs> fragSizeCalcs = Lists.newArrayList();

        for(Map.Entry<String,List<GeneData>> entry : chrGeneMap.entrySet())
        {
            FragmentSizeCalcs fragSizeCalc = new FragmentSizeCalcs(mConfig, mGeneTransCache, mResultsWriter.getFragmentLengthWriter());
            fragSizeCalc.initialise(entry.getKey(), entry.getValue(), requiredFragCount);
            fragSizeCalcs.add(fragSizeCalc);
        }

        final List<Callable> callableList = fragSizeCalcs.stream().collect(Collectors.toList());
        boolean validExecution = TaskExecutor.executeTasks(callableList, mConfig.Threads);

        if(!validExecution)
            return;

        // merge results from all chromosomes
        for(final FragmentSizeCalcs fragSizeCalc : fragSizeCalcs)
        {
            mMaxObservedReadLength = max(mMaxObservedReadLength, fragSizeCalc.getMaxReadLength());
            FragmentSizeCalcs.mergeData(mFragmentLengthDistribution, fragSizeCalc);
        }

        if(mConfig.ApplyFragmentLengthAdjust)
        {
            ISF_LOGGER.info("max observed read length({}) set", mMaxObservedReadLength); // purely for informational purposes
            setConfigFragmentLengthData(mConfig, mFragmentLengthDistribution);
        }

        if (mConfig.WriteFragmentLengths)
        {
            FragmentSizeCalcs.writeFragmentLengths(mConfig, mFragmentLengthDistribution);
        }

        final PerformanceCounter combinedPc = fragSizeCalcs.get(0).getPerformanceCounter();

        for(int i = 1; i < fragSizeCalcs.size(); ++i)
        {
            combinedPc.merge(fragSizeCalcs.get(i).getPerformanceCounter());
        }

        combinedPc.logStats();
    }

    private boolean generateGcRatios(final Map<String,List<GeneData>> chrGeneMap)
    {
        ISF_LOGGER.info("generating GC counts cache");

        // mTranscriptGcRatios.generateExpectedCounts(mChromosome, mGeneDataList);
        final List<GcTranscriptCalculator> taskList = Lists.newArrayList();
        final List<Callable> callableList = Lists.newArrayList();

        mGcTranscriptCalcs.initialiseWriter();

        for(Map.Entry<String,List<GeneData>> entry : chrGeneMap.entrySet())
        {
            GcTranscriptCalculator gcCalcs = new GcTranscriptCalculator(mConfig, mGeneTransCache);
            gcCalcs.initialise(entry.getKey(), entry.getValue(), mGcTranscriptCalcs.getWriter());
            taskList.add(gcCalcs);
            callableList.add(gcCalcs);
        }

        boolean taskStatus = TaskExecutor.executeTasks(callableList, mConfig.Threads);
        mGcTranscriptCalcs.close();
        return taskStatus;
    }

    private boolean generateExpectedCounts(final Map<String,List<GeneData>> chrGeneMap)
    {
        ISF_LOGGER.info("generating expected transcript counts cache");

        final List<ExpressionCacheTask> taskList = Lists.newArrayList();
        final List<Callable> callableList = Lists.newArrayList();

        for(Map.Entry<String,List<GeneData>> entry : chrGeneMap.entrySet())
        {
            ExpressionCacheTask expressionTask = new ExpressionCacheTask(mConfig, mGeneTransCache, mResultsWriter);
            expressionTask.initialise(entry.getKey(), entry.getValue());
            taskList.add(expressionTask);
            callableList.add(expressionTask);
        }

        return TaskExecutor.executeTasks(callableList, mConfig.Threads);
    }

    private boolean countBamReads(final Map<String,List<GeneData>> chrGeneMap)
    {
        ISF_LOGGER.info("basic BAM read counts");

        final List<BamReadCounter> taskList = Lists.newArrayList();
        final List<Callable> callableList = Lists.newArrayList();

        for(Map.Entry<String,List<GeneData>> entry : chrGeneMap.entrySet())
        {
            BamReadCounter bamReaderTask = new BamReadCounter(mConfig);
            bamReaderTask.initialise(entry.getKey(), entry.getValue());
            taskList.add(bamReaderTask);
            callableList.add(bamReaderTask);
        }

        return TaskExecutor.executeTasks(callableList, mConfig.Threads);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final VersionInfo version = new VersionInfo("isofox.version");
        ISF_LOGGER.info("Isofox version: {}", version.version());

        final Options options = createCmdLineOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        if(!validConfigPaths(cmd))
        {
            ISF_LOGGER.error("invalid input files or paths, exiting");
            return;
        }

        IsofoxConfig config = new IsofoxConfig(cmd);

        if(!config.isValid())
        {
            ISF_LOGGER.error("missing config options, exiting");
            return;
        }

        Isofox isofox = new Isofox(config, cmd);
        if(!isofox.runAnalysis())
        {
            ISF_LOGGER.info("Isofox RNA analysis failed");
            return;
        }

        ISF_LOGGER.info("Isofox RNA analysis complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private void logPerformanceStats(final List<PerformanceCounter[]> perfCounters)
    {
        final PerformanceCounter[] combinedPc = perfCounters.get(0);

        for(int i = 1; i < perfCounters.size(); ++i)
        {
            final PerformanceCounter[] chrPCs = perfCounters.get(i);

            for(int j = 0; j < combinedPc.length; ++j)
            {
                combinedPc[j].merge(chrPCs[j]);
            }
        }

        Arrays.stream(combinedPc).forEach(x -> x.logStats());

        if(mConfig.RunPerfChecks)
        {
            // log 10 slowest times and their interval names
            final List<Double> fitTimes = combinedPc[PERF_FIT].getTimes();
            final List<String> fitGenes = combinedPc[PERF_FIT].getTimeNames();

            if(fitTimes.size() >= 10 && fitGenes.size() == fitTimes.size())
            {
                for (int i = fitTimes.size() - 1; i >= fitTimes.size() - 10; --i)
                {
                    ISF_LOGGER.info(String.format("fit times: geneSet(%s) time(%.3f)", fitGenes.get(i), fitTimes.get(i)));
                }
            }
        }
    }

}
