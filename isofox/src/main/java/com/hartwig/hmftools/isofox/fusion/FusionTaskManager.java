package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.HIGH_LOG_COUNT;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.ReadGroup.mergeChimericReadMaps;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.common.utils.TaskExecutor;

public class FusionTaskManager
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final FusionWriter mFusionWriter;
    private final PassingFusions mPassingFusions;

    private final Map<String,List<FusionFragment>> mRealignCandidateMap;
    private final Map<String,Map<String,ReadGroup>> mIncompleteReadGroups; // keyed by chromosome then readId
    private final Map<String,Set<String>> mHardFilteredReadGroups; // keyed by chromosome then readId

    private final PerformanceCounter mPerfCounter;

    public FusionTaskManager(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mPassingFusions = new PassingFusions(config.Fusions.KnownFusions, config.Fusions.CohortFile);

        mRealignCandidateMap = Maps.newHashMap();
        mIncompleteReadGroups = Maps.newHashMap();
        mHardFilteredReadGroups = Maps.newHashMap();

        mPerfCounter = new PerformanceCounter("Fusions");
        mFusionWriter = new FusionWriter(mConfig);
    }

    public FusionFinder createFusionFinder(final String id)
    {
        return new FusionFinder(id, mConfig, mGeneTransCache, mPassingFusions, mFusionWriter);
    }

    public synchronized List<ReadGroup> addIncompleteReadGroup(
            final String chromosome, final Map<String,Map<String,ReadGroup>> chrIncompleteGroups,
            final Map<String,Set<String>> hardFilteredGroups)
    {
        int prevIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();
        int hardFiltered = 0;

        List<ReadGroup> completeGroups = Lists.newArrayList();

        for(Map.Entry<String,Map<String,ReadGroup>> entry : chrIncompleteGroups.entrySet())
        {
            String otherChromosome = entry.getKey();
            Map<String,ReadGroup> newIncompleteGroups = entry.getValue();

            Map<String,ReadGroup> existingGroups = mIncompleteReadGroups.get(otherChromosome);

            if(existingGroups == null)
            {
                Map<String,ReadGroup> chromosomeGroups = mIncompleteReadGroups.get(chromosome);
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

        if(!hardFilteredGroups.isEmpty())
        {
            // use the read-IDs from hard-filtered groups to remove existing incomplete groups
            // otherwise cache them until they can be applied
            for(Map.Entry<String,Set<String>> entry : hardFilteredGroups.entrySet())
            {
                String otherChromosome = entry.getKey();
                Set<String> newHfGroups = entry.getValue();

                Map<String, ReadGroup> existingIncompleteGroups = mIncompleteReadGroups.get(otherChromosome);

                if(existingIncompleteGroups != null)
                {
                    List<String> matched = newHfGroups.stream().filter(x -> existingIncompleteGroups.containsKey(x)).collect(Collectors.toList());
                    matched.forEach(x -> existingIncompleteGroups.remove(x));
                    matched.forEach(x -> newHfGroups.remove(x));
                    hardFiltered += matched.size();
                }

                // store any unmatched HF groups
                Set<String> chrHfGroups = mHardFilteredReadGroups.get(otherChromosome);

                if(chrHfGroups == null)
                {
                    chrHfGroups = Sets.newHashSet();
                    mHardFilteredReadGroups.put(otherChromosome, chrHfGroups);
                }

                chrHfGroups.addAll(newHfGroups);
            }
        }

        int newIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        int partialGroupCount = chrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum();

        ISF_LOGGER.info("chr({}) chimeric groups(partial={} complete={}) total incomplete({} -> {}) hf({})",
                chromosome, partialGroupCount, completeGroups.size(), prevIncomplete, newIncomplete, hardFiltered);

        return completeGroups;
    }

    public synchronized void addRealignCandidateFragments(final Map<String,List<FusionFragment>> racFragments)
    {
        mRealignCandidateMap.putAll(racFragments);
        int totalRacFrags = mRealignCandidateMap.values().stream().mapToInt(x -> x.size()).sum();
        int newRacFrags = racFragments.values().stream().mapToInt(x -> x.size()).sum();

        if(totalRacFrags > HIGH_LOG_COUNT)
        {
            ISF_LOGGER.info("realignable candidate fragments total({}) new({})", totalRacFrags, newRacFrags);
        }
    }


    public synchronized final Map<String,List<FusionFragment>> getRealignCandidateMap() { return mRealignCandidateMap; }

    public void close()
    {
        long unfusedRacFrags = mRealignCandidateMap.values().stream()
                .mapToLong(x -> x.stream().filter(y -> y.assignedFusions().isEmpty()).count()).sum();

        int incompleteGroupCount = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        ISF_LOGGER.info("all fusion tasks complete - unfused RAC frags({}) incompleteGroups({})", unfusedRacFrags, incompleteGroupCount);

        // write any unassigned RAC fragments
        mFusionWriter.writeUnfusedFragments(mRealignCandidateMap);

        if(mConfig.RunPerfChecks)
        {
            List<ReadGroup> incompleteGroups = Lists.newArrayList();

            for(Map.Entry<String,Map<String,ReadGroup>> chrEntry : mIncompleteReadGroups.entrySet())
            {
                String chromosome = chrEntry.getKey();

                if(mConfig.Filters.excludeChromosome(chromosome))
                    continue;

                Map<String, ReadGroup> rgMap = chrEntry.getValue();
                for(ReadGroup readGroup : rgMap.values())
                {
                    if(!mConfig.Filters.SpecificChromosomes.isEmpty())
                    {
                        if(readGroup.Reads.stream().anyMatch(x -> !mConfig.Filters.SpecificChromosomes.contains(x.mateChromosome())))
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

    private boolean skipMissingReads(final List<ReadRecord> reads)
    {
        for(final ReadRecord read : reads)
        {
            if(read.hasSuppAlignment())
            {
                SupplementaryReadData suppData = SupplementaryReadData.from(read.getSuppAlignment());

                if(mConfig.Filters.skipRead(suppData.Chromosome, suppData.Position))
                    return true;

                ISF_LOGGER.debug("read({}) missing supp({})", read, suppData);
            }
        }

        return false;
    }

}
