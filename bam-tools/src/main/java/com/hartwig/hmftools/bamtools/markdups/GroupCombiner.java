package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.chromosomeIndicator;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class GroupCombiner
{
    private final boolean mSinglePartition;
    private final RecordWriter mRecordWriter;

    private final Map<String, PartitionCache> mPartitionCacheMap;
    private final Set<String> mProcessedPartitions;
    private final List<Fragment> mResolvedFragments;

    public GroupCombiner(final RecordWriter recordWriter, final boolean singlePartition)
    {
        mRecordWriter = recordWriter;
        mPartitionCacheMap = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mSinglePartition = singlePartition;
        mResolvedFragments = Lists.newArrayList();
    }

    private PartitionCache getOrCreatePartitionCache(final String chrPartition)
    {
        PartitionCache partitionCache = mPartitionCacheMap.get(chrPartition);

        if(partitionCache == null)
        {
            partitionCache = new PartitionCache(chrPartition);
            mPartitionCacheMap.put(chrPartition, partitionCache);
        }

        return partitionCache;
    }

    public synchronized void processPartitionFragments(
            final String chrPartition, final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList,
            final List<Fragment> supplementaries)
    {
        // cache entries which are expecting further reads, otherwise combine reads from existing fragments and set resolved status

        Map<String, Set<Integer>> impactedPositions = Maps.newHashMap();

        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        for(Fragment fragment : resolvedFragments)
        {
            handleResolvedFragment(fragment, chrPartition, impactedPositions);
        }

        for(Fragment fragment : supplementaries)
        {
            handleSupplementary(fragment, partitionCache);
        }

        for(CandidateDuplicates candidateDuplicates : candidateDuplicatesList)
        {
            handleCandidateDuplicates(candidateDuplicates, partitionCache, impactedPositions);
        }

        // check any candidate duplicate positions & fragments impacted by this set of new fragments
        processImpactedPositions(impactedPositions);

        if(!mSinglePartition)
            partitionComplete(chrPartition);
    }

    private void handleResolvedFragment(final Fragment fragment, final String chrPartition,
            final Map<String, Set<Integer>> impactedPositions)
    {
        if(fragment.allReadsPresent())
            return;

        if(!mSinglePartition && !fragment.hasRemotePartitions())
            return;

        String fragmentChromosome = chromosomeIndicator(fragment.reads().get(0).getReferenceName());

        List<String> remotePartitions = fragment.hasRemotePartitions() ? fragment.remotePartitions() : Lists.newArrayList(chrPartition);

        for(String remotePartition : remotePartitions)
        {
            // local GC stores and applies same chromosome partitions only and vice versa
            if(mSinglePartition && !remotePartition.startsWith(fragmentChromosome)) // only applicable for remote GC
                continue;
            else if(!mSinglePartition && remotePartition.startsWith(fragmentChromosome)) // already applied in local GC
                continue;

            boolean isProcessed = mProcessedPartitions.contains(remotePartition);

            boolean isCurrentPartition = mSinglePartition && remotePartition.equals(chrPartition);

            if(isProcessed || isCurrentPartition)
            {
                setResolvedStatus(remotePartition, fragment, impactedPositions);
            }

            if(!isProcessed || isCurrentPartition)
            {
                // store status against the remote partition, not the current one
                storeResolvedStatus(remotePartition, fragment);
            }
        }
    }

    private void setResolvedStatus(final String chrPartition, final Fragment resolvedFragment,
            final Map<String, Set<Integer>> impactedPositions)
    {
        // when a fragment is resolved and there are unclear fragments for the same read ID, then write them with the resolved status
        PartitionCache partitionCache = mPartitionCacheMap.get(chrPartition);

        if(partitionCache == null)
            return;

        Fragment fragment = partitionCache.IncompleteFragments.get(resolvedFragment.id());

        if(fragment == null)
            return;

        if(fragment.status() == UNCLEAR)
        {
            // remove this position info as well
            // only this fragment, since the others may still turn out to be duplicates
            List<Fragment> fragments = partitionCache.CandidateDuplicatesMap.get(fragment.initialPosition());

            if(fragments != null)
            {
                fragments.remove(fragment);
                addImpactedPosition(impactedPositions, chrPartition, fragment.initialPosition());
            }
        }

        fragment.setStatus(resolvedFragment.status());
        mRecordWriter.writeFragment(fragment);
        partitionCache.IncompleteFragments.remove(resolvedFragment.id());
    }

    private void storeResolvedStatus(final String chrPartition, final Fragment resolvedFragment)
    {
        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);
        partitionCache.FragmentStatus.put(resolvedFragment.id(), resolvedFragment.status());
    }

    public void localSupplementary(final Fragment supplementary, final String chrPartition)
    {
        // not synchronised since only called locally
        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);
        handleSupplementary(supplementary, partitionCache);
    }

    private boolean checkResolvedStatus(final Fragment fragment, final PartitionCache partitionCache)
    {
        FragmentStatus status = partitionCache.FragmentStatus.get(fragment.id());

        if(status != null)
        {
            fragment.setStatus(status);
            mRecordWriter.writeFragment(fragment);
            return true;
        }

        return false;
    }

    private void handleSupplementary(final Fragment supplementary, final PartitionCache partitionCache)
    {
        String remotePartition = !mSinglePartition ? supplementary.remotePartitions().get(0) : partitionCache.ChrPartition;

        boolean processedPartition = mProcessedPartitions.contains(remotePartition);

        boolean storeFragment = !processedPartition || mSinglePartition;

        if(processedPartition || mSinglePartition)
        {
            // check for a resolved status (always stored against the remote partition)
            if(checkResolvedStatus(supplementary, partitionCache))
                return;

            PartitionCache remotePartitionCache = mPartitionCacheMap.get(remotePartition);

            if(remotePartitionCache != null)
            {
                // otherwise check for an unclear fragment to add this to
                Fragment fragment = remotePartitionCache.IncompleteFragments.get(supplementary.id());

                if(fragment != null)
                {
                    fragment.merge(supplementary);
                    return;
                }
            }

            // no read group came through on this remote partition - shouldn't happen, but store the read to be written at the end
            storeFragment = true;
        }

        if(storeFragment)
        {
            // store supplementary for when this partition is processed
            Fragment fragment = partitionCache.IncompleteFragments.get(supplementary.id());

            if(fragment != null)
                fragment.merge(supplementary);
            else
                partitionCache.IncompleteFragments.put(supplementary.id(), supplementary);
        }
    }

    private void handleCandidateDuplicates(
            final CandidateDuplicates candidateDuplicates, final PartitionCache partitionCache, final Map<String,Set<Integer>> impactedPositions)
    {
        // scenarios:
        // - first time this fragment has been seen
        // - existing supplementary, no mate
        // - existing mate, also in an unclear fragment and position group
        // - existing mate with a resolved status

        List<Fragment> resolvedFragments = Lists.newArrayList();
        boolean foundMatches = false;

        // first apply a resolved status if known and remove from position fragments
        int i = 0;
        while(i < candidateDuplicates.Fragments.size())
        {
            Fragment unclearFragment = candidateDuplicates.Fragments.get(i);

            // check for a resolved status (always stored against the remote partition)
            if(checkResolvedStatus(unclearFragment, partitionCache))
            {
                foundMatches = true;
                resolvedFragments.add(unclearFragment);
                candidateDuplicates.Fragments.remove(i);

                // no need ot add to impact positions since this is a new item
            }
            else
            {
                ++i;
            }
        }

        String chrPartition = partitionCache.ChrPartition;

        // next check for existing supplementaries or matched unclear fragments from amongst those previously stored
        for(Fragment unclearFragment : candidateDuplicates.Fragments)
        {
            if(!mSinglePartition && !unclearFragment.hasRemotePartitions())
                continue;

            // for the local GC, only unclear fragment in the same partition are checked for resolving
            List<String> remotePartitions = !mSinglePartition ? Lists.newArrayList(unclearFragment.remotePartitions()) : Lists.newArrayList(chrPartition);

            for(String remotePartition : remotePartitions)
            {
                PartitionCache remotePartitionCache = mPartitionCacheMap.get(remotePartition);

                if(remotePartitionCache == null)
                    continue;

                Fragment matchedFragment = remotePartitionCache.IncompleteFragments.get(unclearFragment.id());

                if(matchedFragment == null)
                    continue;

                // add the cached fragment's reads
                unclearFragment.merge(matchedFragment);

                // remove it from the remote partition and note it needs re-evaluation
                remotePartitionCache.IncompleteFragments.remove(unclearFragment.id());

                // need to remove from its remote position group as well, and then deal with that
                if(matchedFragment.status() == UNCLEAR)
                {
                    foundMatches = true;

                    List<Fragment> posFragments =
                            remotePartitionCache.CandidateDuplicatesMap.get(matchedFragment.initialPosition());

                    if(posFragments != null)
                    {
                        posFragments.remove(matchedFragment);

                        // this depletes the existing position's fragment but leaves them with a status of UNCLEAR, possibly NONE
                        addImpactedPosition(impactedPositions, remotePartition, matchedFragment.initialPosition());
                    }
                }
            }
        }

        if(foundMatches)
        {
            // should be able to resolve some or all of the fragments at this position
            resolveUnclearFragments(candidateDuplicates.Fragments, partitionCache);
        }

        // only store this position's fragment if there are further remote reads expected
        if(!candidateDuplicates.Fragments.isEmpty())
        {
            for(Fragment fragment : candidateDuplicates.Fragments)
            {
                partitionCache.IncompleteFragments.put(fragment.id(), fragment);
            }

            partitionCache.CandidateDuplicatesMap.put(candidateDuplicates.Position, candidateDuplicates.Fragments);
        }
    }

    private void processImpactedPositions(final Map<String, Set<Integer>> impactedPositions)
    {
        for(Map.Entry<String, Set<Integer>> entry : impactedPositions.entrySet())
        {
            PartitionCache partitionCache = mPartitionCacheMap.get(entry.getKey());

            if(partitionCache == null)
                continue;

            for(Integer position : entry.getValue())
            {
                List<Fragment> positionFragments = partitionCache.CandidateDuplicatesMap.get(position);

                if(positionFragments == null)
                    continue;

                resolveUnclearFragments(positionFragments, partitionCache);

                if(positionFragments.isEmpty())
                    partitionCache.CandidateDuplicatesMap.remove(position);
            }
        }
    }

    private void resolveUnclearFragments(final List<Fragment> positionFragments, final PartitionCache partitionCache)
    {
        if(positionFragments.isEmpty())
            return;

        List<Fragment> resolvedFragments = Lists.newArrayList();

        positionFragments.forEach(x -> x.setStatus(UNSET)); // reset before evaluation

        classifyFragments(positionFragments, resolvedFragments, null);

        for(Fragment fragment : resolvedFragments)
        {
            mRecordWriter.writeFragment(fragment);
            mResolvedFragments.add(fragment);

            if(partitionCache != null)
                partitionCache.IncompleteFragments.remove(fragment.id());

            // store resolved status if expecting more reads
            handleResolvedFragment(fragment, partitionCache.ChrPartition, null);

            // clear from other unclear position fragments
        }
    }

    private static void addImpactedPosition(final Map<String, Set<Integer>> impactedPositions, final String chrPartition, int position)
    {
        if(impactedPositions == null)
            return;

        Set<Integer> positions = impactedPositions.get(chrPartition);

        if(positions == null)
            impactedPositions.put(chrPartition, Sets.newHashSet(position));
        else
            positions.add(position);
    }

    private void partitionComplete(final String chrPartition)
    {
        mProcessedPartitions.add(chrPartition);

        // write unmatched local supplementaries
        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        // remove fragment status since no longer required
        partitionCache.FragmentStatus.clear();

        if(!mSinglePartition && (mProcessedPartitions.size() % 100) == 0)
        {
            int cachedFragments = mPartitionCacheMap.values().stream().mapToInt(x -> x.incompleteFragments()).sum();
            int cachedResolved = mPartitionCacheMap.values().stream().mapToInt(x -> x.resolvedFragments()).sum();

            BM_LOGGER.info("group cache partitions processed({}) cached fragments({}) resolved({})",
                    mProcessedPartitions.size(), cachedFragments, cachedResolved);
        }
    }

    public void gatherPartitionFragments(
            final String chrPartition, final BaseRegion partitionRegion,
            final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList)
    {
        if(!mSinglePartition)
        {
            BM_LOGGER.error("remote GC called local function for partition({})", chrPartition);
            return;
        }

        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        for(Fragment fragment : mResolvedFragments)
        {
            fragment.setRemotePartitions(partitionRegion);

            if(fragment.hasRemotePartitions())
                resolvedFragments.add(fragment);
        }

        for(Map.Entry<Integer, List<Fragment>> entry : partitionCache.CandidateDuplicatesMap.entrySet())
        {
            List<Fragment> unclearFragments = entry.getValue();
            unclearFragments.forEach(x -> x.setRemotePartitions(partitionRegion));
            candidateDuplicatesList.add(new CandidateDuplicates(entry.getKey(), unclearFragments));
        }

        for(Fragment fragment : partitionCache.IncompleteFragments.values())
        {
            if(fragment.status() == SUPPLEMENTARY)
            {
                mRecordWriter.writeFragment(fragment);
                BM_LOGGER.debug("supplementary({}) in local GC unmatched", fragment);
            }
        }

        partitionCache.clear();
        mProcessedPartitions.add(chrPartition);
        mResolvedFragments.clear();
    }

    public void handleRemaining()
    {
        int cachedFragments = mPartitionCacheMap.values().stream().mapToInt(x -> x.incompleteFragments()).sum();
        int cachedResolved = mPartitionCacheMap.values().stream().mapToInt(x -> x.resolvedFragments()).sum();

        BM_LOGGER.info("final group cache partitions({}) incomplete({}) resolved({})",
                mProcessedPartitions.size(), cachedFragments, cachedResolved);

        for(PartitionCache partitionCache : mPartitionCacheMap.values())
        {
            for(Fragment fragment : partitionCache.IncompleteFragments.values())
            {
                mRecordWriter.writeFragment(fragment);
            }
        }

        mPartitionCacheMap.clear();
    }

    public void reset()
    {
        mPartitionCacheMap.clear();
        mProcessedPartitions.clear();
    }

    public PartitionCache getPartitionCache(final String chrPartition) { return mPartitionCacheMap.get(chrPartition); }
}