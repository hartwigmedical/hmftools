package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.HIGH_LOG_COUNT;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.GeneCollection;

public class ChromosomeFusions
{
    private final String mChromosome;
    private final IsofoxConfig mConfig;
    private final FusionTaskManager mFusionTaskManager;
    private final ChimericReadTracker mChimericReadTracker;

    private final FusionFinder mFusionFinder;
    private final ChimericStats mChimericStats;
    private final PerformanceCounter mPerfCounter;

    public ChromosomeFusions(
            final IsofoxConfig config, final String chromosome,
            final FusionTaskManager fusionManager, final ChimericReadTracker chimericReadTracker, final PerformanceCounter perfCounter)
    {
        mConfig = config;
        mChromosome = chromosome;
        mFusionTaskManager = fusionManager;
        mChimericReadTracker = chimericReadTracker;
        mPerfCounter = perfCounter;
        
        mChimericStats = new ChimericStats();

        mFusionFinder = mFusionTaskManager.createFusionFinder(mChromosome);

        mChimericReadTracker.setKnownSpliteSites(fusionManager.hardFilteredCache().getKnownSpliteSites());

    }

    public ChimericStats chimericStats() { return mChimericStats; }
    
    public void onGeneCollectionComplete(final GeneCollection geneCollection, final BaseDepth baseDepth)
    {
        if(!mConfig.runFunction(FUSIONS) || mFusionFinder == null)
            return;

        if(!mPerfCounter.isRunning())
            mPerfCounter.start();
        else
            mPerfCounter.resume();

        // pass any complete chimeric read groups to the fusion finder
        // and add to this any groups which are now complete (ie which were partially complete before)
        // cache any incomplete groups, either for later gene collections or from other chromosomes
        final List<FusionReadGroup> completeReadGroups = mFusionFinder.processNewChimericReadGroups(
                geneCollection, baseDepth,mChimericReadTracker.getReadMap());

        mChimericStats.merge(mChimericReadTracker.getStats());

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
                mChromosome, geneCollection.id(), mChimericReadTracker.extractJunctionRacFragments());

        mFusionFinder.processLocalReadGroups(completeReadGroups);

        if(highCount)
        {
            mFusionFinder.clearState(false);
            System.gc();
        }

        mPerfCounter.pause();

    }
    
    public void onChromosomeComplete()
    {
        if(mChimericStats.ChimericJunctions > HIGH_LOG_COUNT)
        {
            ISF_LOGGER.info("chr({}) chimeric data: {}", mChromosome, mChimericStats);
        }

        mPerfCounter.stop();

        mPerfCounter.start();

        // handle fragments spanning multiple chromosomes

        Map<String, Set<String>> chrHardFilteredIds = mChimericReadTracker.getHardFilteredReadIds();

        // organise incomplete reads into the chromosomes which they link to
        final Map<String,Map<String,FusionReadGroup>> chrIncompleteReadsGroups = mFusionFinder.extractIncompleteReadGroups(
                mChromosome, chrHardFilteredIds);

        final List<FusionReadGroup> interChromosomalGroups = mFusionTaskManager.addIncompleteReadGroup(
                mChromosome, chrIncompleteReadsGroups, chrHardFilteredIds);

        if(!interChromosomalGroups.isEmpty())
        {
            mFusionFinder.processInterChromosomalReadGroups(interChromosomalGroups);
        }

        mPerfCounter.stop();

        mChimericReadTracker.clearAll();

        mFusionFinder.logPerfCounters();
        mFusionFinder.clearState(true);

    }
}
