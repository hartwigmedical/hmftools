package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionReadGroup.mergeChimericReadMaps;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class FusionTaskManager
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final FusionWriter mFusionWriter;
    private final PassingFusions mPassingFusions;

    private final RacFragmentCache mRacFragmentCache;
    private final Map<String,Map<String, FusionReadGroup>> mIncompleteReadGroups; // keyed by chromosome then readId
    private final Map<String,Set<String>> mHardFilteredReadGroups; // keyed by chromosome then readId

    public FusionTaskManager(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mPassingFusions = new PassingFusions(config.Fusions.KnownFusions, config.Fusions.CohortFile);

        mRacFragmentCache = new RacFragmentCache();
        mIncompleteReadGroups = Maps.newHashMap();
        mHardFilteredReadGroups = Maps.newHashMap();

        mFusionWriter = new FusionWriter(mConfig);
    }

    public FusionFinder createFusionFinder(final String id)
    {
        return new FusionFinder(id, mConfig, mGeneTransCache, mRacFragmentCache, mPassingFusions, mFusionWriter);
    }

    public final RacFragmentCache racFragmentCache() { return mRacFragmentCache; }

    public synchronized List<FusionReadGroup> addIncompleteReadGroup(
            final String chromosome, final Map<String,Map<String, FusionReadGroup>> chrIncompleteGroups)
    {
        int prevIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();
        int hardFiltered = 0;

        List<FusionReadGroup> completeGroups = Lists.newArrayList();

        for(Map.Entry<String,Map<String, FusionReadGroup>> entry : chrIncompleteGroups.entrySet())
        {
            String otherChromosome = entry.getKey();
            Map<String, FusionReadGroup> newIncompleteGroups = entry.getValue();

            Map<String, FusionReadGroup> existingGroups = mIncompleteReadGroups.get(otherChromosome);

            if(existingGroups == null)
            {
                Map<String, FusionReadGroup> chromosomeGroups = mIncompleteReadGroups.get(chromosome);
                if(chromosomeGroups == null)
                {
                    chromosomeGroups = Maps.newHashMap();
                    mIncompleteReadGroups.put(chromosome, chromosomeGroups);
                }

                chromosomeGroups.putAll(newIncompleteGroups);

                ISF_LOGGER.debug("added chromosomes({} & {}) newGroups({})",
                        chromosome, otherChromosome, newIncompleteGroups.size());
            }
            else
            {
                Set<String> chrHfGroups = mHardFilteredReadGroups.get(otherChromosome);

                if(chrHfGroups != null)
                {
                    List<String> matched = chrHfGroups.stream().filter(x -> newIncompleteGroups.containsKey(x)).collect(Collectors.toList());
                    matched.forEach(x -> newIncompleteGroups.remove(x));
                    matched.forEach(x -> chrHfGroups.remove(x));
                    hardFiltered += matched.size();
                }

                mergeChimericReadMaps(existingGroups, completeGroups, newIncompleteGroups);

                ISF_LOGGER.debug("combined chromosomes({} & {}) existing({}) new({}) complete({})",
                        chromosome, otherChromosome, existingGroups.size(), newIncompleteGroups.size(), completeGroups.size());
            }
        }

        int newIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        int partialGroupCount = chrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum();

        ISF_LOGGER.info("chr({}) chimeric groups(partial={} complete={}) total incomplete({} -> {}) hf({})",
                chromosome, partialGroupCount, completeGroups.size(), prevIncomplete, newIncomplete, hardFiltered);

        return completeGroups;
    }

    public RacFragmentCache getRacFragmentCache() { return mRacFragmentCache; }

    public synchronized void addRacFragments(final String chromosome, int geneCollectionId, final JunctionRacFragments racFragments)
    {
        mRacFragmentCache.addRacFragments(chromosome, geneCollectionId, racFragments);
    }

    public void close()
    {
        int incompleteGroupCount = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        ISF_LOGGER.info("all fusion tasks complete - RAC frags({} assigned={} groups={}) incompleteGroups({})",
                mRacFragmentCache.totalFragmentCount(), mRacFragmentCache.assignedFragmentCount(),
                mRacFragmentCache.totalGroupCount(), incompleteGroupCount);

        // write any unassigned RAC fragments
        if(mConfig.Fusions.WriteChimericReads || mConfig.Fusions.WriteChimericFragments)
            mFusionWriter.writeUnfusedFragments(mRacFragmentCache.getUnassignedFragments());

        if(mConfig.RunPerfChecks)
        {
            List<FusionReadGroup> incompleteGroups = Lists.newArrayList();

            for(Map.Entry<String,Map<String, FusionReadGroup>> chrEntry : mIncompleteReadGroups.entrySet())
            {
                String chromosome = chrEntry.getKey();

                if(mConfig.Filters.excludeChromosome(chromosome))
                    continue;

                Map<String, FusionReadGroup> rgMap = chrEntry.getValue();
                for(FusionReadGroup readGroup : rgMap.values())
                {
                    if(!mConfig.Filters.SpecificChromosomes.isEmpty())
                    {
                        if(readGroup.Reads.stream().anyMatch(x -> !mConfig.Filters.SpecificChromosomes.contains(x.MateChromosome)))
                            continue;
                    }

                    if(!skipMissingReads(readGroup.Reads))
                    {
                        incompleteGroups.add(readGroup);
                    }
                }
            }

            mFusionWriter.writeIncompleteGroupReads(incompleteGroups);
        }

        mFusionWriter.close();
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
