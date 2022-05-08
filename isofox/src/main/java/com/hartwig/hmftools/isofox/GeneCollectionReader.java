package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.RETAINED_INTRONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.GeneReadData.createGeneReadData;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_FIT;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_FUSIONS;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_GC_ADJUST;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_NOVEL_LOCATIONS;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_READS;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.PERF_TOTAL;
import static com.hartwig.hmftools.isofox.common.PerformanceTracking.logMemory;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.isofox.common.CommonUtils.getChromosomeLength;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.HIGH_LOG_COUNT;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.common.FragmentType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.PerformanceTracking;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.expression.ExpectedCountsCache;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesData;
import com.hartwig.hmftools.isofox.expression.TranscriptExpression;
import com.hartwig.hmftools.isofox.expression.GeneCollectionSummary;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.adjusts.GcTranscriptCalculator;
import com.hartwig.hmftools.isofox.fusion.ChimericStats;
import com.hartwig.hmftools.isofox.fusion.FusionFinder;
import com.hartwig.hmftools.isofox.fusion.FusionTaskManager;
import com.hartwig.hmftools.isofox.fusion.FusionReadGroup;
import com.hartwig.hmftools.isofox.results.GeneResult;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.isofox.results.TranscriptResult;

public class GeneCollectionReader implements Callable
{
    private final String mChromosome;
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;
    private final ResultsWriter mResultsWriter;

    private final BamFragmentAllocator mBamFragmentAllocator;
    private final TranscriptExpression mExpTransRates;
    private final GcTranscriptCalculator mTranscriptGcRatios;
    private final ExpectedCountsCache mExpectedCountsCache;

    private final List<GeneData> mGeneDataList;
    private int mCollectionId;
    private int mCurrentGeneIndex;
    private int mGenesProcessed;

    // fusion state cached across all gene collections
    private final FusionTaskManager mFusionTaskManager;
    private final FusionFinder mFusionFinder;
    private final ChimericStats mChimericStats;

    // cache of results
    private final List<GeneCollectionSummary> mGeneCollectionSummaryData;
    private int mEnrichedGenesFragmentCount;
    private final int[] mCombinedFragmentCounts;
    private final GcRatioCounts mNonEnrichedGcRatioCounts;
    private int mTotalReadsProcessed;
    private final GcRatioCounts mGcRatioCounts;

    private TaskType mCurrentTaskType;
    private boolean mIsValid;

    private final PerformanceCounter[] mPerfCounters;

    public GeneCollectionReader(
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

        mBamFragmentAllocator = new BamFragmentAllocator(mConfig, resultsWriter);
        mBamFragmentAllocator.registerKnownFusionPairs(mGeneTransCache);

        mGcRatioCounts = mBamFragmentAllocator.getGcRatioCounts();
        mExpTransRates = mConfig.ExpCountsFile != null ? new TranscriptExpression(mConfig, mExpectedCountsCache, resultsWriter) : null;
        mTranscriptGcRatios = transcriptGcCalcs;

        mGeneCollectionSummaryData = Lists.newArrayList();
        mEnrichedGenesFragmentCount = 0;
        mTotalReadsProcessed = 0;
        mCombinedFragmentCounts = new int[typeAsInt(FragmentType.MAX)];
        mNonEnrichedGcRatioCounts = new GcRatioCounts();
        mChimericStats = new ChimericStats();

        mFusionTaskManager = fusionManager;
        mFusionFinder = mFusionTaskManager != null ? mFusionTaskManager.createFusionFinder(mChromosome) : null;

        mPerfCounters = PerformanceTracking.createPerfCounters();

        mIsValid = true;
    }

    public String chromosome() { return mChromosome; }
    public final List<GeneCollectionSummary> getGeneCollectionSummaryData() { return mGeneCollectionSummaryData; }
    public final GcRatioCounts getGcRatioCounts() { return mGcRatioCounts; }

    public final ChimericStats getChimericStats() { return mChimericStats; }
    public boolean isValid() { return mIsValid; }
    public int totalReadCount() { return mTotalReadsProcessed; }

    public void setTaskType(TaskType taskType) { mCurrentTaskType = taskType; }

    @Override
    public Long call()
    {
        if(mCurrentTaskType == null)
        {
            ISF_LOGGER.error(" no chromosome-gene task set for execution");
            return (long)0;
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

        return (long)1; // return value not used
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

        boolean genesFiltered = !mConfig.Filters.RestrictedGeneIds.isEmpty() || !mConfig.Filters.SpecificRegions.isEmpty();

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
                    int endOfChromosome = (int)getChromosomeLength(mChromosome, mConfig.RefGenVersion);
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

                if(mConfig.runFunction(FUSIONS))
                    ISF_LOGGER.debug("chr({}) chimeric data: {}", mChromosome, mChimericStats);
            }
        }

        if(mConfig.runFunction(FUSIONS))
        {
            if(nextLogCount > 100)
            {
                ISF_LOGGER.info("chr({}) chimeric data: {}", mChromosome, mChimericStats);
            }

            mPerfCounters[PERF_FUSIONS].stop();

            mPerfCounters[PERF_FUSIONS].start();

            // handle fragments spanning multiple chromosomes

            Map<String,Set<String>> chrHardFilteredIds = mBamFragmentAllocator.getChimericReadTracker().getHardFilteredReadIds();

            // organise incomplete reads into the chromosomes which they link to
            final Map<String,Map<String,FusionReadGroup>> chrIncompleteReadsGroups = mFusionFinder.extractIncompleteReadGroups(
                    mChromosome, chrHardFilteredIds);

            final List<FusionReadGroup> interChromosomalGroups = mFusionTaskManager.addIncompleteReadGroup(
                    mChromosome, chrIncompleteReadsGroups, chrHardFilteredIds);

            if(!interChromosomalGroups.isEmpty())
            {
                mFusionFinder.processInterChromosomalReadGroups(interChromosomalGroups);
            }

            mPerfCounters[PERF_FUSIONS].stop();

            mBamFragmentAllocator.getChimericReadTracker().clearAll();

            mFusionFinder.logPerfCounters();
            mFusionFinder.clearState(true);
        }

        if(mGeneDataList.size() > 10)
        {
            ISF_LOGGER.info("chr({}) processing complete", mChromosome);
        }

        logMemory(mConfig, String.format("chr(%s)-Complete", mChromosome));
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
            for(int i = 0; i < mCombinedFragmentCounts.length; ++i)
            {
                mCombinedFragmentCounts[i] += geneCollection.getCounts()[i];
            }

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

            final Map<Integer,String> transIdMap = Maps.newHashMap();
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

        if(!mConfig.Filters.EnrichedGeneIds.isEmpty())
        {
            int enrichedGeneFragments = geneCollection.genes().stream()
                    .anyMatch(x -> mConfig.Filters.EnrichedGeneIds.contains(x.GeneData.GeneId))
                    ? geneCollection.getCounts()[typeAsInt(TOTAL)] : 0;

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

        for(int i = 0; i < mCombinedFragmentCounts.length; ++i)
        {
            mCombinedFragmentCounts[i] += geneCollection.getCounts()[i];
        }

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
        if(!mConfig.runFunction(FUSIONS) || mFusionFinder == null)
            return;

        if(!mPerfCounters[PERF_FUSIONS].isRunning())
            mPerfCounters[PERF_FUSIONS].start();
        else
            mPerfCounters[PERF_FUSIONS].resume();

        // pass any complete chimeric read groups to the fusion finder
        // and add to this any groups which are now complete (ie which were partially complete before)
        // cache any incomplete groups, either for later gene collections or from other chromosomes
        final List<FusionReadGroup> completeReadGroups = mFusionFinder.processNewChimericReadGroups(
                geneCollection, mBamFragmentAllocator.getBaseDepth(),
                mBamFragmentAllocator.getChimericReadTracker().getReadMap());

        mChimericStats.merge(mBamFragmentAllocator.getChimericReadTracker().getStats());

        boolean highCount = completeReadGroups.size() >= HIGH_LOG_COUNT;
        if(highCount)
        {
            int nonSuppGroups = (int)completeReadGroups.stream().filter(x -> x.size() == 2).count();
            ISF_LOGGER.info("chr({}) genes({}) region({} - {}) found {} local chimeric read groups (non-supp={}), stats({})",
                    mChromosome, geneCollection.geneNames(),
                    geneCollection.getNonGenicPositions()[SE_START], geneCollection.getNonGenicPositions()[SE_END],
                    completeReadGroups.size(), nonSuppGroups, mChimericStats);
        }

        mFusionTaskManager.addRacFragments(
                mChromosome, geneCollection.id(), mBamFragmentAllocator.getChimericReadTracker().extractJunctionRacFragments());

        mFusionFinder.processLocalReadGroups(completeReadGroups);

        if(highCount)
        {
            mFusionFinder.clearState(false);
            System.gc();
        }

        mPerfCounters[PERF_FUSIONS].pause();
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

            final Map<Integer,String> transIdMap = Maps.newHashMap();
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

    public int getEnrichedGenesFragmentCount() { return mEnrichedGenesFragmentCount; }
    public int[] getCombinedCounts() { return mCombinedFragmentCounts; }
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
                        .filter(x -> x.GeneData.GeneId.equals(transResult.Trans.GeneId))
                        .map(x -> x.GeneData)
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
