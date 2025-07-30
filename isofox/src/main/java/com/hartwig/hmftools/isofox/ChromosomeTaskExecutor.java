package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.RETAINED_INTRONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.common.CommonUtils.getChromosomeLength;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.GeneReadData.createGeneReadData;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_FIT;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_FUSIONS;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_GC_ADJUST;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_NOVEL_LOCATIONS;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_READS;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_TOTAL;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findUniqueBases;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.adjusts.GcTranscriptCalculator;
import com.hartwig.hmftools.isofox.common.FragmentTypeCounts;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.PerformanceTracking;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.expression.ExpectedCountsCache;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesData;
import com.hartwig.hmftools.isofox.expression.GeneCollectionSummary;
import com.hartwig.hmftools.isofox.expression.TranscriptExpression;
import com.hartwig.hmftools.isofox.fusion.ChimericStats;
import com.hartwig.hmftools.isofox.fusion.ChromosomeFusions;
import com.hartwig.hmftools.isofox.fusion.FusionTaskManager;
import com.hartwig.hmftools.isofox.results.GeneResult;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.isofox.results.TranscriptResult;

public class ChromosomeTaskExecutor implements Callable<Void>
{
    private final String mChromosome;
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;
    private final ResultsWriter mResultsWriter;

    private final FragmentAllocator mBamFragmentAllocator;
    private final TranscriptExpression mExpTransRates;
    private final GcTranscriptCalculator mTranscriptGcRatios;
    private final ExpectedCountsCache mExpectedCountsCache;

    private final List<GeneData> mGeneDataList;
    private int mCollectionId;
    private int mCurrentGeneIndex;
    private int mGenesProcessed;

    // fusion state cached across all gene collections
    private final ChromosomeFusions mChromosomeFusions;

    // cache of results
    private final List<GeneCollectionSummary> mGeneCollectionSummaryData;
    private long mEnrichedGenesFragmentCount;
    private final FragmentTypeCounts mCombinedFragmentCounts;
    private final GcRatioCounts mNonEnrichedGcRatioCounts;
    private int mTotalReadsProcessed;
    private final GcRatioCounts mGcRatioCounts;

    private TaskType mCurrentTaskType;
    private final boolean mIsValid;

    private final PerformanceCounter[] mPerfCounters;

    public ChromosomeTaskExecutor(
            final IsofoxConfig config, final String chromosome, final List<GeneData> geneDataList,
            final EnsemblDataCache geneTransCache, final ResultsWriter resultsWriter, final FusionTaskManager fusionManager,
            final ExpectedCountsCache expectedCountsCache, final GcTranscriptCalculator transcriptGcCalcs)
    {
        mConfig = config;
        mChromosome = chromosome;
        mGeneTransCache = geneTransCache;
        mResultsWriter = resultsWriter;

        mGeneDataList = geneDataList;
        mCollectionId = 0;

        mCurrentGeneIndex = 0;
        mCurrentTaskType = null;

        mExpectedCountsCache = expectedCountsCache;

        mBamFragmentAllocator = new FragmentAllocator(mConfig, resultsWriter);
        mBamFragmentAllocator.registerKnownFusionPairs(mGeneTransCache);

        mGcRatioCounts = mBamFragmentAllocator.getGcRatioCounts();
        mExpTransRates = mConfig.ExpCountsFile != null ? new TranscriptExpression(mConfig, mExpectedCountsCache, resultsWriter) : null;
        mTranscriptGcRatios = transcriptGcCalcs;

        mGeneCollectionSummaryData = Lists.newArrayList();
        mEnrichedGenesFragmentCount = 0;
        mTotalReadsProcessed = 0;
        mCombinedFragmentCounts = new FragmentTypeCounts();
        mNonEnrichedGcRatioCounts = new GcRatioCounts();

        mPerfCounters = PerformanceTracking.createPerfCounters();

        mChromosomeFusions = mConfig.runFunction(FUSIONS) ? new ChromosomeFusions(
                config, chromosome, fusionManager, mBamFragmentAllocator.getChimericReadTracker(), mPerfCounters[PERF_FUSIONS]) : null;

        mIsValid = true;
    }

    public String chromosome() { return mChromosome; }
    public final List<GeneCollectionSummary> getGeneCollectionSummaryData() { return mGeneCollectionSummaryData; }
    public final GcRatioCounts getGcRatioCounts() { return mGcRatioCounts; }

    public final ChimericStats getChimericStats() { return mChromosomeFusions.chimericStats(); }
    public boolean isValid() { return mIsValid; }
    public int totalReadCount() { return mTotalReadsProcessed; }

    public void setTaskType(TaskType taskType) { mCurrentTaskType = taskType; }

    @Override
    public Void call()
    {
        if(mCurrentTaskType == null)
        {
            ISF_LOGGER.error(" no chromosome-gene task set for execution");
            return null;
        }

        switch(mCurrentTaskType)
        {
            case TRANSCRIPT_COUNTS:
                assignTranscriptCounts();
                break;

            case APPLY_GC_ADJUSTMENT:
                applyGcAdjustToTranscriptAllocations();
                break;

            default:
                break;
        }

        return null;
    }

    public void assignTranscriptCounts()
    {
        if(mGeneDataList.size() > 10)
        {
            ISF_LOGGER.info("chr({}) processing {} genes", mChromosome, mGeneDataList.size());
        }

        mCurrentGeneIndex = 0;
        final List<GeneData> overlappingGenes = Lists.newArrayList();
        int nextLogCount = 100;
        int lastGeneCollectionEndPosition = 1;

        boolean genesFiltered = !mConfig.Filters.RestrictedGeneIds.isEmpty() || mConfig.Filters.SpecificChrRegions.hasFilters();

        while(mCurrentGeneIndex < mGeneDataList.size())
        {
            mCurrentGeneIndex = findNextOverlappingGenes(mGeneDataList, mCurrentGeneIndex, overlappingGenes);

            final List<GeneReadData> geneReadDataList = createGeneReadData(overlappingGenes, mGeneTransCache);

            GeneCollection geneCollection = new GeneCollection(mCollectionId++, geneReadDataList);
            geneCollection.markEnrichedAndExcludedGenes(mConfig, mGeneTransCache);

            if(!genesFiltered) // reads will be taken from the previous gene collection's end
            {
                geneCollection.setNonGenicPosition(SE_START, lastGeneCollectionEndPosition);

                if(mCurrentGeneIndex < mGeneDataList.size())
                {
                    final GeneData nextGeneData = mGeneDataList.get(mCurrentGeneIndex);
                    geneCollection.setNonGenicPosition(SE_END, nextGeneData.GeneStart - 1);
                }
                else
                {
                    int endOfChromosome = (int) getChromosomeLength(mChromosome, mConfig.RefGenVersion);
                    int endNonGenicPosition = max(geneCollection.getNonGenicPositions()[SE_START] + 1, endOfChromosome - 1000);
                    geneCollection.setNonGenicPosition(SE_END, endNonGenicPosition);
                    geneCollection.setEndOfChromosome();
                }
            }
            else
            {
                // the buffer is to be able to test out pre and post gene region reads
                if(lastGeneCollectionEndPosition == 1)
                {
                    geneCollection.setNonGenicPosition(SE_START, geneCollection.regionBounds()[SE_START] - 10000);
                }
                else
                {
                    geneCollection.setNonGenicPosition(SE_START, lastGeneCollectionEndPosition);
                }

                if(mCurrentGeneIndex < mGeneDataList.size())
                {
                    final GeneData nextGeneData = mGeneDataList.get(mCurrentGeneIndex);
                    geneCollection.setNonGenicPosition(SE_END, nextGeneData.GeneStart - 1);
                }
                else
                {
                    geneCollection.setNonGenicPosition(SE_END, geneCollection.regionBounds()[SE_END] + 10000);
                }
            }

            mPerfCounters[PERF_TOTAL].start();

            analyseBamReads(geneCollection);

            mPerfCounters[PERF_TOTAL].stop();

            ISF_LOGGER.debug("chr({}) gene({}) processed({} of {})",
                    mChromosome, geneCollection.geneNames(10), mCurrentGeneIndex, mGeneDataList.size());

            mGenesProcessed += geneCollection.genes().size();
            mTotalReadsProcessed = mBamFragmentAllocator.totalReadCount();

            lastGeneCollectionEndPosition = geneCollection.regionBounds()[SE_END] + 1;

            if(mGenesProcessed >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.info("chr({}) processed {} of {} genes", mChromosome, mGenesProcessed, mGeneDataList.size());
            }
        }

        if(mChromosomeFusions != null)
            mChromosomeFusions.onChromosomeComplete();

        if(mGeneDataList.size() > 10)
        {
            ISF_LOGGER.info("chr({}) processing complete", mChromosome);
        }
    }

    public static int findNextOverlappingGenes(
            final List<GeneData> geneDataList, int currentIndex, final List<GeneData> overlappingGenes)
    {
        overlappingGenes.clear();

        while(currentIndex < geneDataList.size())
        {
            GeneData geneData = geneDataList.get(currentIndex);

            if(overlappingGenes.isEmpty()
            || overlappingGenes.stream().anyMatch(x -> positionsOverlap(geneData.GeneStart, geneData.GeneEnd, x.GeneStart, x.GeneEnd)))
            {
                overlappingGenes.add(geneData);
                ++currentIndex;
            }
            else
            {
                break;
            }
        }

        return currentIndex;
    }

    private void analyseBamReads(final GeneCollection geneCollection)
    {
        // cache reference bases for comparison with read bases
        if(mConfig.RefGenomeFile != null)
        {
            for(RegionReadData region : geneCollection.getExonRegions())
            {
                final String regionRefBases = mConfig.RefGenome.getBaseString(region.chromosome(), region.start(), region.end());
                region.setRefBases(regionRefBases);
            }

            findUniqueBases(geneCollection.getExonRegions());
        }

        // start the read region at the previous gene collection's end if known
        int[] geneRegionPositions;

        if(mConfig.runFunction(FUSIONS))
        {
            geneRegionPositions = geneCollection.getNonGenicPositions();
        }
        else
        {
            geneRegionPositions = new int[] { geneCollection.regionBounds()[SE_START] - 100, geneCollection.regionBounds()[SE_END] + 100 };
        }

        if(geneRegionPositions[SE_START] >= geneRegionPositions[SE_END])
        {
            ISF_LOGGER.warn("invalid geneCollection({}) region({} -> {})",
                    geneCollection.geneNames(), geneRegionPositions[SE_START], geneRegionPositions[SE_END]);
            return;
        }

        final ChrBaseRegion geneRegion = new ChrBaseRegion(geneCollection.chromosome(), geneRegionPositions);

        mPerfCounters[PERF_READS].start();
        mBamFragmentAllocator.produceBamCounts(geneCollection, geneRegion);
        mPerfCounters[PERF_READS].stop();

        postBamReadTranscriptCounts(geneCollection);
        postBamReadNovelLocations(geneCollection);
        postBamReadFusions(geneCollection);

        mBamFragmentAllocator.clearCache(); // free up resources for this gene collection
    }

    private void postBamReadTranscriptCounts(final GeneCollection geneCollection)
    {
        if(mConfig.runStatisticsOnly())
        {
            mCombinedFragmentCounts.combine(geneCollection.fragmentTypeCounts());
            return;
        }

        if(!mConfig.runFunction(TRANSCRIPT_COUNTS))
            return;

        GeneCollectionSummary geneCollectionSummary = new GeneCollectionSummary(
                geneCollection.chrId(), geneCollection.geneIds(), geneCollection.geneNames(), mBamFragmentAllocator.getTransComboData());

        mGeneCollectionSummaryData.add(geneCollectionSummary);

        if(ISF_LOGGER.isDebugEnabled())
        {
            double allCategoryTotals = mBamFragmentAllocator.getTransComboData().stream()
                    .mapToDouble(x -> x.fragmentCount()).sum();

            if(allCategoryTotals > 0)
            {
                double transCategoryTotals = mBamFragmentAllocator.getTransComboData().stream()
                        .filter(x -> !x.transcriptIds().isEmpty())
                        .mapToDouble(x -> x.fragmentCount()).sum();

                ISF_LOGGER.debug(String.format("genes(%s) catCounts(all=%.2f trans=%.1f)",
                        geneCollection.geneNames(), allCategoryTotals, transCategoryTotals));
            }

            ISF_LOGGER.debug("chr({}) gene({}) transCombo(gene={} total={})",
                    mChromosome, geneCollection.geneNames(10), mBamFragmentAllocator.getTransComboData().size(),
                    mGeneCollectionSummaryData.stream().mapToInt(x -> x.TransCategoryCounts.size()).sum());
        }

        if(mExpTransRates != null)
        {
            ExpectedRatesData expRatesData = null;

            mPerfCounters[PERF_FIT].start();

            final Map<Integer, String> transIdMap = Maps.newHashMap();
            geneCollection.getTranscripts().forEach(x -> transIdMap.put(x.TransId, x.TransName));
            mExpTransRates.runTranscriptEstimation(transIdMap, geneCollectionSummary, expRatesData, false);

            mPerfCounters[PERF_FIT].stop();
        }

        for(GeneReadData geneReadData : geneCollection.genes())
        {
            collectResults(geneCollection, geneCollectionSummary, geneReadData);

            if(mConfig.WriteExonData)
            {
                geneReadData.getTranscripts().forEach(x -> mResultsWriter.writeExonData(geneReadData, x));
            }
        }

        if(mConfig.WriteSpliceJunctions)
        {
            mResultsWriter.writeSpliceJunctionData(geneCollection);
        }

        if(!mConfig.Filters.EnrichedGeneIds.isEmpty())
        {
            long enrichedGeneFragments = geneCollection.genes().stream()
                    .anyMatch(x -> mConfig.Filters.EnrichedGeneIds.contains(x.GeneData.GeneId))
                    ? geneCollection.fragmentTypeCounts().typeCount(TOTAL) : 0;

            if(enrichedGeneFragments > 0)
            {
                mEnrichedGenesFragmentCount += enrichedGeneFragments;
            }
            else
            {
                if(mBamFragmentAllocator.getGeneGcRatioCounts() != null)
                    mNonEnrichedGcRatioCounts.mergeRatioCounts(mBamFragmentAllocator.getGeneGcRatioCounts().getCounts());
            }
        }
        else
        {
            // take them all
            if(mBamFragmentAllocator.getGeneGcRatioCounts() != null)
                mNonEnrichedGcRatioCounts.mergeRatioCounts(mBamFragmentAllocator.getGeneGcRatioCounts().getCounts());
        }

        mCombinedFragmentCounts.combine(geneCollection.fragmentTypeCounts());

        geneCollectionSummary.allocateResidualsToGenes();
        mResultsWriter.writeGeneCollectionData(geneCollection);

        if(!mConfig.applyGcBiasAdjust())
            geneCollectionSummary.TransCategoryCounts.clear();
    }

    private void postBamReadNovelLocations(final GeneCollection geneCollection)
    {
        if(!mConfig.runFunction(ALT_SPLICE_JUNCTIONS) && !mConfig.runFunction(RETAINED_INTRONS))
            return;

        mPerfCounters[PERF_NOVEL_LOCATIONS].start();
        mBamFragmentAllocator.annotateNovelLocations();

        if(mConfig.WriteSpliceSiteData)
        {
            mBamFragmentAllocator.getSpliceSiteCounter().writeSpliceSiteData(geneCollection);
        }

        mPerfCounters[PERF_NOVEL_LOCATIONS].stop();
    }

    private void postBamReadFusions(final GeneCollection geneCollection)
    {
        if(mChromosomeFusions == null)
            return;

        mChromosomeFusions.onGeneCollectionComplete(geneCollection, mBamFragmentAllocator.getBaseDepth());
    }

    public void applyGcAdjustment()
    {
        mPerfCounters[PERF_GC_ADJUST].start();
        mTranscriptGcRatios.generateGcCountsFromFit(mGeneCollectionSummaryData);
        mPerfCounters[PERF_GC_ADJUST].pause();
    }

    private void applyGcAdjustToTranscriptAllocations()
    {
        mPerfCounters[PERF_GC_ADJUST].resume();

        for(final GeneCollectionSummary geneSummaryData : mGeneCollectionSummaryData)
        {
            final double[] gcAdjustments = mTranscriptGcRatios.getGcRatioAdjustments();
            geneSummaryData.applyGcAdjustments(gcAdjustments);

            final Map<Integer, String> transIdMap = Maps.newHashMap();
            geneSummaryData.TranscriptResults.forEach(x -> transIdMap.put(x.Trans.TransId, x.Trans.TransName));
            mExpTransRates.runTranscriptEstimation(transIdMap, geneSummaryData, null, true);

            geneSummaryData.setFitAllocations();
            geneSummaryData.allocateResidualsToGenes();
        }

        mPerfCounters[PERF_GC_ADJUST].stop();
    }

    private void collectResults(
            final GeneCollection geneCollection, final GeneCollectionSummary geneCollectionSummary, final GeneReadData geneReadData)
    {
        for(final TranscriptData transData : geneReadData.getTranscripts())
        {
            final TranscriptResult results = new TranscriptResult(geneCollection, geneReadData, transData, mConfig.FragmentSizeData);

            geneCollectionSummary.TranscriptResults.add(results);
        }

        GeneResult geneResult = new GeneResult(geneCollection, geneReadData);
        geneCollectionSummary.GeneResults.add(geneResult);

        geneCollectionSummary.setFitAllocations();
        geneCollectionSummary.assignLowMapQualityFragments();
        geneCollectionSummary.TranscriptResults.forEach(x -> x.setPreGcFitAllocation(x.getFitAllocation()));
    }

    public long getEnrichedGenesFragmentCount() { return mEnrichedGenesFragmentCount; }
    public FragmentTypeCounts getCombinedCounts() { return mCombinedFragmentCounts; }
    public GcRatioCounts getNonEnrichedGcRatioCounts() { return mNonEnrichedGcRatioCounts; }

    public void writeResults()
    {
        for(final GeneCollectionSummary geneCollectionResult : mGeneCollectionSummaryData)
        {
            for(final GeneResult geneResult : geneCollectionResult.GeneResults)
            {
                mResultsWriter.writeGeneResult(geneResult);
            }

            for(final TranscriptResult transResult : geneCollectionResult.TranscriptResults)
            {
                final GeneData geneData = geneCollectionResult.GeneResults.stream()
                        .filter(x -> x.Gene.GeneId.equals(transResult.Trans.GeneId))
                        .map(x -> x.Gene)
                        .findFirst().orElse(null);

                mResultsWriter.writeTranscriptResults(geneData, transResult);
            }
        }
    }

    public PerformanceCounter[] getPerfCounters()
    {
        return mPerfCounters;
    }
}
