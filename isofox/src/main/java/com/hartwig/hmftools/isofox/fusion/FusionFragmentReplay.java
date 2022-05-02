package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class FusionFragmentReplay
{
    private final IsofoxConfig mConfig;
    private final List<FusionFinder> mFusionTasks;

    private final EnsemblDataCache mGeneTransCache;
    private final RacFragmentCache mRacFragmentCache;

    private final FusionWriter mFusionWriter;
    private final PassingFusions mPassingFusions;

    public FusionFragmentReplay(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mPassingFusions = new PassingFusions(config.Fusions.KnownFusions, config.Fusions.CohortFile);
        mFusionWriter = new FusionWriter(mConfig);
        mRacFragmentCache = new RacFragmentCache();

        mFusionTasks = Lists.newArrayList();
    }

    public FusionFinder createFusionFinder(final String id)
    {
        return new FusionFinder(id, mConfig, mGeneTransCache, mRacFragmentCache, mPassingFusions, mFusionWriter);
    }

    public void run()
    {
        List<FusionReadGroup> cachedReadGroups = ChimericReadCache.loadChimericReads(mConfig.Fusions.ChimericReadsFile);

        processCachedFragments(cachedReadGroups);
    }

    private static final int LOG_COUNT = 100000;

    public void processCachedFragments(final List<FusionReadGroup> readGroups)
    {
        // convert any set of valid reads into a fragment, and then process these in groups by chromosomal pair
        ISF_LOGGER.info("processing {} chimeric read groups", readGroups.size());

        int invalidFragments = 0;
        int skipped = 0;
        int fragments = 0;
        int missingSuppReads = 0;

        int readGroupCount = 0;

        final Map<String,List<FusionFragment>> chrPairFragments = Maps.newHashMap();

        for(FusionReadGroup readGroup : readGroups)
        {
            ++readGroupCount;

            if(readGroupCount > 0 && (readGroupCount % LOG_COUNT) == 0)
            {
                ISF_LOGGER.info("processed {} chimeric read groups", readGroupCount);
            }

            boolean isComplete = readGroup.isComplete();
            final List<FusionRead> reads = readGroup.Reads;

            if(reads.stream().anyMatch(x -> mConfig.Filters.skipRead(x.MateChromosome, x.MatePosStart)))
            {
                ++skipped;
                continue;
            }

            if(!isComplete)
            {
                if(readGroup.hasSuppAlignment())
                {
                    if(skipMissingReads(reads))
                    {
                        ++skipped;
                        continue;
                    }

                    ++missingSuppReads;
                }

                ++invalidFragments;
            }
            else
            {
                FusionFragment fragment = new FusionFragment(readGroup);

                if(fragment.type() == FusionFragmentType.UNKNOWN)
                {
                    ++invalidFragments;
                    continue;
                }

                ++fragments;

                final String chrPair = formChromosomePair(fragment.chromosomes()[SE_START], fragment.chromosomes()[SE_END]);
                List<FusionFragment> fragmentList = chrPairFragments.get(chrPair);

                if(fragmentList == null)
                {
                    chrPairFragments.put(chrPair, Lists.newArrayList(fragment));
                }
                else
                {
                    fragmentList.add(fragment);
                }
            }
        }

        // allocate fusion pairs evenly amongst threads (if multi-thread)
        mFusionTasks.clear();

        for(int taskId = 0; taskId < max(mConfig.Threads, 1); ++taskId)
        {
            mFusionTasks.add(createFusionFinder(String.valueOf(taskId)));
        }

        for(List<FusionFragment> chrPairFrags : chrPairFragments.values())
        {
            // allocate the next chr-pair fragment batch to the task with the least
            FusionFinder leastAllocated = null;
            int minAllocated = 0;

            for(FusionFinder fusionTask : mFusionTasks)
            {
                if(minAllocated == 0 || fusionTask.getFragments().size() < minAllocated)
                {
                    leastAllocated = fusionTask;
                    minAllocated = fusionTask.getFragments().size();
                    if(minAllocated == 0)
                        break;
                }
            }

            leastAllocated.getFragments().addAll(chrPairFrags);
        }

        ISF_LOGGER.info("chimeric groups({} skipped={} invalid={} miss={} candidates={}) chrPairs({}) tasks({})",
                readGroups.size(), skipped, invalidFragments, missingSuppReads, fragments,
                chrPairFragments.size(), mFusionTasks.size());

        if(mFusionTasks.isEmpty())
        {
            ISF_LOGGER.warn("no fusion tasks created");
            return;
        }
        else
        {
            final List<Callable> callableList = mFusionTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
            logPerformanceStats();
        }

        mFusionWriter.close();

        ISF_LOGGER.info("fusion calling complete");
    }

    private void logPerformanceStats()
    {
        if(!mConfig.Fusions.PerformanceStats)
            return;

        if(!ISF_LOGGER.isDebugEnabled() && mConfig.Functions.size() > 1)
            return;

        final PerformanceCounter perfCounter = mFusionTasks.get(0).getPerfCounter();

        for (int i = 1; i < mFusionTasks.size(); ++i)
        {
            perfCounter.merge(mFusionTasks.get(i).getPerfCounter());
        }

        perfCounter.logStats();
    }

    private boolean skipMissingReads(final List<FusionRead> reads)
    {
        for(final FusionRead read : reads)
        {
            if(read.HasSuppAlignment)
            {
                if(mConfig.Filters.skipRead(read.SuppData.Chromosome, read.SuppData.Position))
                    return true;

                ISF_LOGGER.debug("read({}) missing supp({})", read, read.SuppData);
            }
        }

        return false;
    }

}
