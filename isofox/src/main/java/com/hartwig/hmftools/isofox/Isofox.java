package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.rna.RnaStatistics.LOW_COVERAGE_THRESHOLD;
import static com.hartwig.hmftools.common.rna.RnaStatistics.SPLICE_GENE_THRESHOLD;
import static com.hartwig.hmftools.common.sigs.SigUtils.convertToPercentages;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.APP_NAME;
import static com.hartwig.hmftools.isofox.IsofoxConstants.PANEL_LOW_COVERAGE_FACTOR;
import static com.hartwig.hmftools.isofox.IsofoxConstants.PRIORITISED_CHROMOSOMES;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.NEO_EPITOPES;
import static com.hartwig.hmftools.isofox.IsofoxFunction.READ_COUNTS;
import static com.hartwig.hmftools.isofox.TaskType.APPLY_GC_ADJUSTMENT;
import static com.hartwig.hmftools.isofox.TaskType.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs.setConfigFragmentLengthData;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.writeReadGcRatioCounts;
import static com.hartwig.hmftools.isofox.expression.TranscriptExpression.calcTpmFactors;
import static com.hartwig.hmftools.isofox.expression.TranscriptExpression.setTranscriptsPerMillion;
import static com.hartwig.hmftools.isofox.results.SummaryStats.createSummaryStats;

import java.io.IOException;
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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.adjusts.GcTranscriptCalculator;
import com.hartwig.hmftools.isofox.common.BamReadCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.isofox.common.FragmentTypeCounts;
import com.hartwig.hmftools.isofox.common.PerformanceTracking;
import com.hartwig.hmftools.isofox.cram.CramAndValidateCommands;
import com.hartwig.hmftools.isofox.cram.StageRunner;
import com.hartwig.hmftools.isofox.expression.ExpectedCountsCache;
import com.hartwig.hmftools.isofox.expression.GeneCollectionSummary;
import com.hartwig.hmftools.isofox.expression.PanelTpmNormaliser;
import com.hartwig.hmftools.isofox.fusion.ChimericStats;
import com.hartwig.hmftools.isofox.fusion.FusionTaskManager;
import com.hartwig.hmftools.isofox.neo.NeoEpitopeReader;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

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
    private final PerformanceTracking mPerfTracking;

    public Isofox(final IsofoxConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        mResultsWriter = new ResultsWriter(mConfig);
        mPerfTracking = new PerformanceTracking(mConfig);

        mGeneTransCache = new EnsemblDataCache(configBuilder);

        if(!mConfig.Filters.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.setRestrictedGeneIdList(mConfig.Filters.RestrictedGeneIds);
        }

        mGeneTransCache.setRequiredData(true, false, false, mConfig.CanonicalTranscriptOnly);
        mGeneTransCache.load(false);

        mConfig.Filters.buildGeneRegions(mGeneTransCache);

        mExpectedCountsCache = mConfig.ExpCountsFile != null || mConfig.applyGcBiasAdjust() ? new ExpectedCountsCache(mConfig) : null;

        mGcTranscriptCalcs = mConfig.applyGcBiasAdjust() ? new GcTranscriptCalculator(mConfig) : null;

        mFusionTaskManager = mConfig.runFunction(FUSIONS) ? new FusionTaskManager(mConfig, mGeneTransCache) : null;

        mMaxObservedReadLength = 0;
        mFragmentLengthDistribution = Lists.newArrayList();
    }

    public boolean runAnalysis()
    {
        long startTimeMs = System.currentTimeMillis();

        // all other routines split work by chromosome
        Map<String,List<GeneData>> chrGeneMap = getChromosomeGeneLists();

        if(chrGeneMap.isEmpty())
        {
            ISF_LOGGER.error("no chromosome tasks created");
            return false;
        }

        if(mConfig.runFunction(NEO_EPITOPES))
        {
            NeoEpitopeReader neReader = new NeoEpitopeReader(mConfig, mGeneTransCache);
            neReader.calcFragmentSupport();
            return true;
        }

        if(mConfig.runFunction(READ_COUNTS))
        {
            boolean status = countBamReads(chrGeneMap);
            mResultsWriter.close();
            return status;
        }

        // create CRAM file
        try {
            ISF_LOGGER.info("Cram export");
            ISF_LOGGER.info("BAM file:{}", mConfig.BamFile);
            ISF_LOGGER.info("Output dir for CRAM:{}", mConfig.OutputDir);
            new StageRunner().run(mConfig.BamFile, mConfig.OutputDir);
        } catch (IOException e) {
            System.exit(1);
            throw new RuntimeException(e);

        }

        // BAM processing for the key routines - novel junctions, fusions and gene expression
        if(!allocateBamFragments(chrGeneMap))
            return false;

        ISF_LOGGER.info("Isofox complete, mins({})", runTimeMinsStr(startTimeMs));
        System.exit(1);
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

        final List<ChromosomeTaskExecutor> chrTasks = Lists.newArrayList();
        final List<Callable> callableList = Lists.newArrayList();
        final List<String> chromosomes = Lists.newArrayList();

        // process any enriched genes first, then add the rest in order of decreasing length
        PRIORITISED_CHROMOSOMES.forEach(x -> chromosomes.add(mConfig.RefGenVersion.versionedChromosome(x)));

        Arrays.stream(HumanChromosome.values())
                .map(chromosome -> mConfig.RefGenVersion.versionedChromosome(chromosome.toString()))
                .filter(chromosome -> !chromosomes.contains(chromosome))
                .forEach(chromosome -> chromosomes.add(chromosome));

        for(String chromosome : chromosomes)
        {
            List<GeneData> geneDataList = chrGeneMap.get(chromosome);

            if(geneDataList == null)
                continue;

            ChromosomeTaskExecutor bamReaderTask = new ChromosomeTaskExecutor(
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
            // post-processing for summary stats and gene expression data
            processBamFragments(chrTasks, callableList);
        }

        if(mConfig.runFunction(FUSIONS))
        {
            // extract all chimeric reads and associated data for fusion calling
            ChimericStats chimericStats = new ChimericStats();
            chrTasks.forEach(x -> chimericStats.merge(x.getChimericStats()));
            ISF_LOGGER.info("overall chimeric stats: {} inv={}", chimericStats, chimericStats.Inversions);

            mFusionTaskManager.close();
        }

        final List<PerformanceCounter[]> perfCounters = chrTasks.stream().map(x -> x.getPerfCounters()).collect(Collectors.toList());
        chrTasks.clear();

        mPerfTracking.logPerformanceStats(perfCounters);
        return true;
    }

    private void processBamFragments(final List<ChromosomeTaskExecutor> chrTasks, final List<Callable> callableList)
    {
        FragmentTypeCounts totalFragmentCounts = new FragmentTypeCounts();

        long enrichedGeneFragCount = 0;
        GcRatioCounts nonEnrichedGcRatioCounts = new GcRatioCounts();

        for(ChromosomeTaskExecutor chrTask : chrTasks)
        {
            totalFragmentCounts.combine(chrTask.getCombinedCounts());

            enrichedGeneFragCount += chrTask.getEnrichedGenesFragmentCount();

            nonEnrichedGcRatioCounts.mergeRatioCounts(chrTask.getNonEnrichedGcRatioCounts().getCounts());
        }

        if(mConfig.applyGcBiasAdjust())
        {
            applyGcAdjustments(chrTasks, callableList, nonEnrichedGcRatioCounts);
        }

        if(mConfig.WriteGcData)
        {
            GcRatioCounts combinedGcRatioCounts = new GcRatioCounts();
            chrTasks.forEach(x -> combinedGcRatioCounts.mergeRatioCounts(x.getGcRatioCounts().getCounts()));

            writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "ALL", combinedGcRatioCounts.getCounts(), false);
            double[] percentData = new double[combinedGcRatioCounts.size()];

            copyVector(combinedGcRatioCounts.getCounts(), percentData);
            convertToPercentages(percentData);
            writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "ALL_PERC", percentData, true);

            if(!mConfig.Filters.EnrichedGeneIds.isEmpty())
                writeReadGcRatioCounts(mResultsWriter.getReadGcRatioWriter(), "NON_ENRICHED", nonEnrichedGcRatioCounts.getCounts(), false);

            if(mConfig.applyGcBiasAdjust())
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

        double[] tpmFactors = calcTpmFactors(geneSummaryData, mConfig.Filters.EnrichedGeneIds);

        PanelTpmNormaliser panelTpmNormaliser = new PanelTpmNormaliser(mConfig.PanelTpmNormFile);

        int spliceGeneCount = 0;

        for(ChromosomeTaskExecutor chrTask : chrTasks)
        {
            setTranscriptsPerMillion(chrTask.getGeneCollectionSummaryData(), tpmFactors);

            panelTpmNormaliser.applyNormalisation(chrTask.getGeneCollectionSummaryData());

            spliceGeneCount += chrTask.getGeneCollectionSummaryData().stream().mapToInt(x -> x.spliceGenesCount()).sum();

            chrTask.writeResults();
        }

        // write summary statistics
        if(mConfig.runFunction(IsofoxFunction.TRANSCRIPT_COUNTS) || mConfig.runFunction(IsofoxFunction.STATISTICS))
        {
            double medianGCRatio = nonEnrichedGcRatioCounts.getPercentileRatio(0.5);

            int lowCoverageThreshold = LOW_COVERAGE_THRESHOLD;
            int splicedGeneThreshold = SPLICE_GENE_THRESHOLD;

            if(!mConfig.Filters.RestrictedGeneIds.isEmpty() && !mConfig.PanelTpmNormFile.isEmpty())
            {
                // could be adjusted for the specific panel
                double panelGeneCoverage = mConfig.Filters.RestrictedGeneIds.size() / 37000.0; // total gene count
                lowCoverageThreshold = (int)(lowCoverageThreshold * panelGeneCoverage * PANEL_LOW_COVERAGE_FACTOR);
                splicedGeneThreshold = (int)(panelGeneCoverage * SPLICE_GENE_THRESHOLD);
            }

            final RnaStatistics summaryStats = createSummaryStats(
                    totalFragmentCounts, enrichedGeneFragCount, spliceGeneCount,
                    medianGCRatio, mFragmentLengthDistribution, mMaxObservedReadLength > 0 ? mMaxObservedReadLength : mConfig.ReadLength,
                    lowCoverageThreshold, splicedGeneThreshold);

            mResultsWriter.writeSummaryStats(summaryStats);
        }

        mResultsWriter.close();
    }

    private void applyGcAdjustments(
            final List<ChromosomeTaskExecutor> chrTasks, final List<Callable> callableList, final GcRatioCounts actualGcCounts)
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
        if(!mConfig.Filters.SpecificChrRegions.hasFilters())
            return mGeneTransCache.getChrGeneDataMap();

        final Map<String,List<GeneData>> chrGeneMap = Maps.newHashMap();

        if(mConfig.Filters.SpecificChrRegions.Regions.isEmpty())
        {
            mGeneTransCache.getChrGeneDataMap().entrySet().stream()
                    .filter(x -> mConfig.Filters.SpecificChrRegions.includeChromosome(x.getKey()))
                    .filter(x -> !x.getValue().isEmpty())
                    .forEach(x -> chrGeneMap.put(x.getKey(), x.getValue()));
        }
        else
        {
            for(ChrBaseRegion region : mConfig.Filters.SpecificChrRegions.Regions)
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

        if(mConfig.WriteFragmentLengths)
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

    private boolean countBamReads(final Map<String,List<GeneData>> chrGeneMap)
    {
        ISF_LOGGER.info("basic BAM read counts");

        final List<BamReadCounter> taskList = Lists.newArrayList();
        final List<Callable> callableList = Lists.newArrayList();

        for(Map.Entry<String,List<GeneData>> entry : chrGeneMap.entrySet())
        {
            BamReadCounter bamReaderTask = new BamReadCounter(mConfig, mResultsWriter);
            bamReaderTask.initialise(entry.getKey(), entry.getValue());
            taskList.add(bamReaderTask);
            callableList.add(bamReaderTask);
        }

        return TaskExecutor.executeTasks(callableList, mConfig.Threads);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        IsofoxConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        IsofoxConfig config = new IsofoxConfig(configBuilder);

        if(!config.isValid())
        {
            ISF_LOGGER.error("missing config options, exiting");
            System.exit(1);
        }

        Isofox isofox = new Isofox(config, configBuilder);

        if(!isofox.runAnalysis())
        {
            ISF_LOGGER.info("Isofox RNA analysis failed");
            System.exit(1);
        }
    }
}
