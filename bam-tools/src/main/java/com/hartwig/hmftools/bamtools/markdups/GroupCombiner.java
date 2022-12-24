package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class GroupCombiner
{
    private final boolean mSinglePartition;
    private final RecordWriter mRecordWriter;

    private final Map<String,PartitionCache> mPartitionCacheMap;
    private final Set<String> mProcessedPartitions;

    private class PartitionCache
    {
        // fragment status from resolved fragments, keyed by chromosome-partition then readId
        public final Map<String,FragmentStatus> FragmentStatus;

        // incomplete fragments (unclear or supplmentaries), keyed by chromosome-partition then readId
        public final Map<String,Fragment> IncompleteFragments;

        // positions with candidate duplicate fragments, keyed by chromosome-partition then initial fragment coordinate position
        public final Map<Integer,List<Fragment>> IncompleteFragmentPositions;

        public PartitionCache()
        {
            FragmentStatus = Maps.newHashMap();
            IncompleteFragments = Maps.newHashMap();
            IncompleteFragmentPositions = Maps.newHashMap();
        }

        public int incompleteFragments() { return IncompleteFragments.size(); }
        public int resolvedFragments() { return FragmentStatus.size(); }
    }

    public GroupCombiner(final RecordWriter recordWriter, final boolean singlePartition)
    {
        mRecordWriter = recordWriter;
        mPartitionCacheMap = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mSinglePartition = singlePartition;
    }

    public synchronized void processPartitionFragments(
            final String chrPartition, final List<Fragment> resolvedFragments, final List<PositionFragments> incompletePositionFragments,
            final List<Fragment> supplementaries)
    {
        // cache entries which are expecting further reads, otherwise combine reads from existing fragments and set resolved status

        Map<String,Set<Integer>> impactedPositions = Maps.newHashMap();

        for(Fragment fragment : resolvedFragments)
        {
            handleResolvedFragment(fragment, chrPartition, impactedPositions);
        }

        for(Fragment fragment : supplementaries)
        {
            handleSupplementary(fragment, chrPartition);
        }

        for(PositionFragments positionFragments : incompletePositionFragments)
        {
            handleIncompleteFragments(positionFragments, chrPartition, impactedPositions);
        }

        // process any cached position fragments from this set of new fragments
        processImpactedPositions(impactedPositions);

        if(!mSinglePartition)
            partitionComplete(chrPartition);
    }

    private void handleResolvedFragment(final Fragment fragment, final String chrPartition, final Map<String,Set<Integer>> impactedPositions)
    {
        if(!mSinglePartition)
        {
            for(String remotePartition : fragment.remotePartitions())
            {
                if(mProcessedPartitions.contains(remotePartition))
                {
                    setResolvedStatus(remotePartition, fragment, impactedPositions);
                }
                else
                {
                    // store status for when this partition is processed
                    storeResolvedStatus(chrPartition, fragment);
                }
            }
        }
        else
        {
            setResolvedStatus(chrPartition, fragment, impactedPositions);
            storeResolvedStatus(chrPartition, fragment);
        }
    }

    private void setResolvedStatus(final String chrPartition, final Fragment resolvedFragment, final Map<String,Set<Integer>> impactedPositions)
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
            List<Fragment> fragments = partitionCache.IncompleteFragmentPositions.get(fragment.initialPosition());

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
        PartitionCache partitionCache = mPartitionCacheMap.get(chrPartition);

        if(partitionCache == null)
        {
            partitionCache = new PartitionCache();
            mPartitionCacheMap.put(chrPartition, partitionCache);
        }

        partitionCache.FragmentStatus.put(resolvedFragment.id(), resolvedFragment.status());
    }

    public void processSupplementary(final Fragment supplementary, final String chrPartition)
    {
        // not synchronised since only called locally
        handleSupplementary(supplementary, chrPartition);
    }

    private void handleSupplementary(final Fragment supplementary, final String chrPartition)
    {
        String remotePartition = !mSinglePartition ? supplementary.remotePartitions().get(0) : chrPartition;

        boolean processedPartition = mProcessedPartitions.contains(remotePartition);

        boolean storeFragment = !processedPartition || mSinglePartition;

        if(processedPartition || mSinglePartition)
        {
            // check for a resolved status
            PartitionCache partitionCache = mPartitionCacheMap.get(remotePartition);

            if(partitionCache != null)
            {
                FragmentStatus status = partitionCache.FragmentStatus.get(supplementary.id());

                if(status != null)
                {
                    supplementary.setStatus(status);
                    mRecordWriter.writeFragment(supplementary);
                    return;
                }

                // otherwise check for an unclear fragment to add this to
                Fragment fragment = partitionCache.IncompleteFragments.get(supplementary.id());

                if(fragment != null)
                {
                    supplementary.reads().forEach(x -> fragment.addRead(x));
                    return;
                }
            }

            // no read group came through on this remote partition - shouldn't happen, but store the read to be written at the end
            storeFragment = true;
        }

        if(storeFragment)
        {
            // store supplementary for when this partition is processed
            PartitionCache partitionCache = mPartitionCacheMap.get(chrPartition);

            if(partitionCache == null)
            {
                partitionCache = new PartitionCache();
                mPartitionCacheMap.put(chrPartition, partitionCache);
            }

            Fragment fragment = partitionCache.IncompleteFragments.get(supplementary.id());

            if(fragment != null)
                supplementary.reads().forEach(x -> fragment.addRead(x));
            else
                partitionCache.IncompleteFragments.put(supplementary.id(), supplementary);
        }
    }

    private void handleIncompleteFragments(
            final PositionFragments positionFragments, final String chrPartition, final Map<String,Set<Integer>> impactedPositions)
    {
        // scenarios:
        // - first time this fragment has been seen
        // - existing supplementary, no mate
        // - existing mate, also in an unclear fragment and position group
        // - existing mate with a resolved status

        // steps:
        // - first apply a resolved status if known and remove from position fragments
        // - gather any existing supplementaries and matching unclear fragments
        // - if mate is not present, store for now
        // - cache if awaiting mate reads otherwise combine reads from matched fragment and resolve
        // - add mate reads which will then resolve a set of fragments
        List<Fragment> resolvedFragments = Lists.newArrayList();
        boolean foundMatches = false;

        for(Fragment unclearFragment : positionFragments.Fragments)
        {
            if(!mSinglePartition && !unclearFragment.hasRemotePartitions())
                continue;

            List<String> remotePartitions = !mSinglePartition ? unclearFragment.remotePartitions() : Lists.newArrayList(chrPartition);

            for(String remotePartition : remotePartitions)
            {
                if(mSinglePartition || mProcessedPartitions.contains(chrPartition))
                {
                    // check if this fragment has a resolved status
                    PartitionCache partitionCache = mPartitionCacheMap.get(remotePartition);

                    if(partitionCache == null)
                        continue;

                    FragmentStatus status = partitionCache.FragmentStatus.get(unclearFragment.id());

                    if(status != null)
                    {
                        foundMatches = true;

                        resolvedFragments.add(unclearFragment);
                        unclearFragment.setStatus(status);

                        // can't remove this entry since it may be needed for other reads (eg supplementaries)
                    }
                }
            }
        }

        // process up these resolved fragments
        for(Fragment fragment : resolvedFragments)
        {
            positionFragments.Fragments.remove(fragment);
            mRecordWriter.writeFragment(fragment);

            // no need ot add to impact positions since this is a new item
            // addImpactedPosition(impactedPositions, chrPartition, fragment.initialPosition());
        }

        // next check for existing supplementaries or matched unclear fragments from amongst those previously stored
        for(Fragment unclearFragment : positionFragments.Fragments)
        {
            if(!mSinglePartition && !unclearFragment.hasRemotePartitions())
                continue;

            List<String> remotePartitions = !mSinglePartition ? unclearFragment.remotePartitions() : Lists.newArrayList(chrPartition);

            for(String remotePartition : remotePartitions)
            {
                PartitionCache partitionCache = mPartitionCacheMap.get(remotePartition);

                if(partitionCache == null)
                    continue;

                Fragment matchedFragment = partitionCache.IncompleteFragments.get(unclearFragment.id());

                if(matchedFragment != null)
                {
                    if(matchedFragment.status() == UNCLEAR)
                        foundMatches = true;

                    // add the cached fragment's reads
                    matchedFragment.reads().forEach(x -> unclearFragment.addRead(x));

                    // remove it from the remote partition and note it needs re-evaluation
                    partitionCache.IncompleteFragments.remove(unclearFragment.id());

                    // need to remove from its remote position group as well, and then deal with that
                    List<Fragment> posFragments = partitionCache.IncompleteFragmentPositions.get(matchedFragment.initialPosition());

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
            resolveUnclearFragments(positionFragments.Fragments, null);
        }

        // only store this position's fragment if there are further remote reads expected
        if(!positionFragments.Fragments.isEmpty())
        {
            PartitionCache partitionCache = mPartitionCacheMap.get(chrPartition);

            if(partitionCache == null)
            {
                partitionCache = new PartitionCache();
                mPartitionCacheMap.put(chrPartition, partitionCache);
            }

            for(Fragment fragment : positionFragments.Fragments)
            {
                partitionCache.IncompleteFragments.put(fragment.id(), fragment);
            }

            partitionCache.IncompleteFragmentPositions.put(positionFragments.Position, positionFragments.Fragments);
        }
    }

    private void processImpactedPositions(final Map<String,Set<Integer>> impactedPositions)
    {
        for(Map.Entry<String,Set<Integer>> entry : impactedPositions.entrySet())
        {
            PartitionCache partitionCache = mPartitionCacheMap.get(entry.getKey());

            if(partitionCache == null)
                continue;

            for(Integer position : entry.getValue())
            {
                List<Fragment> positionFragments = partitionCache.IncompleteFragmentPositions.get(position);

                if(positionFragments == null)
                    continue;

                resolveUnclearFragments(positionFragments, partitionCache);

                if(positionFragments.isEmpty())
                    partitionCache.IncompleteFragmentPositions.remove(position);
            }
        }
    }

    private void resolveUnclearFragments(final List<Fragment> positionFragments, final PartitionCache partitionCache)
    {
        List<Fragment> resolvedFragments = Lists.newArrayList();

        positionFragments.forEach(x -> x.setStatus(UNSET)); // reset before evaluation

        classifyFragments(positionFragments, resolvedFragments, null);

        for(Fragment fragment : resolvedFragments)
        {
            // positionFragments.remove(fragment);
            mRecordWriter.writeFragment(fragment);

            if(partitionCache != null)
                partitionCache.IncompleteFragments.remove(fragment.id());
        }
    }

    private static void addImpactedPosition(final Map<String,Set<Integer>> impactedPositions, final String chrPartition, int position)
    {
        Set<Integer> positions = impactedPositions.get(chrPartition);

        if(positions == null)
            impactedPositions.put(chrPartition, Sets.newHashSet(position));
        else
            positions.add(position);
    }

    private void partitionComplete(final String chrPartition)
    {
        mProcessedPartitions.add(chrPartition);

        if((mProcessedPartitions.size() % 100) == 0)
        {
            int cachedFragments = mPartitionCacheMap.values().stream().mapToInt(x -> x.incompleteFragments()).sum();
            int cachedResolved = mPartitionCacheMap.values().stream().mapToInt(x -> x.resolvedFragments()).sum();

            BM_LOGGER.info("group cache partitions processed({}) cached fragments({}) resolved({})",
                    mProcessedPartitions.size(), cachedFragments, cachedResolved);
        }
    }

    public Map<Integer,List<Fragment>> incompleteFragmentPositions(final String chrParitition)
    {
        PartitionCache partitionCache = mPartitionCacheMap.get(chrParitition);
        return partitionCache != null ? partitionCache.IncompleteFragmentPositions : Collections.emptyMap();
    }

    public List<Fragment> incompleteFragments(final String chrParitition)
    {
        PartitionCache partitionCache = mPartitionCacheMap.get(chrParitition);
        return partitionCache != null ? partitionCache.IncompleteFragments.values().stream().collect(Collectors.toList()) : Collections.emptyList();
    }

    public List<Fragment> unmatchedSupplementaries(final String chrParitition)
    {
        PartitionCache partitionCache = mPartitionCacheMap.get(chrParitition);
        return partitionCache != null ? partitionCache.IncompleteFragments.values().stream()
                .filter(x -> x.status() == SUPPLEMENTARY)
                .collect(Collectors.toList()) : Collections.emptyList();
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
                // fragment.setStatus(NONE); // leave their status as-is to assist with evaluation
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
}
