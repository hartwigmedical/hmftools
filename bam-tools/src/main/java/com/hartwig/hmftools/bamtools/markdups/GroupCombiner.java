package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class GroupCombiner
{
    private final RecordWriter mRecordWriter;

    // fragment status from resolved fragments, keyed by chromosome-partition then readId
    private final Map<String,Map<String,FragmentStatus>> mFragmentStatus;

    // incomplete fragments (unclear or supplmentaries), keyed by chromosome-partition then readId
    private final Map<String,Map<String,Fragment>> mIncompleteFragments;

    // positions with candidate duplicate fragments, keyed by chromosome-partition then initial fragment coordinate position
    private final Map<String,Map<Integer,List<Fragment>>> mIncompleteFragmentPositions;

    private final Set<String> mProcessedPartitions;
    private final boolean mSinglePartition;

    public GroupCombiner(final RecordWriter recordWriter, final boolean singlePartition)
    {
        mRecordWriter = recordWriter;
        mFragmentStatus = Maps.newHashMap();
        mIncompleteFragments = Maps.newHashMap();
        mIncompleteFragmentPositions = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mSinglePartition = singlePartition;
    }

    public Map<String,Map<String,FragmentStatus>> fragmentStatusMap() { return mFragmentStatus; }
    public Map<String,Map<String,Fragment>> incompleteFragmentsMap() { return mIncompleteFragments; }
    public Map<String,Map<Integer,List<Fragment>>> incompleteFragmentPositionsMap() { return mIncompleteFragmentPositions; }

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
        Map<String,Fragment> readIdMap = mIncompleteFragments.get(chrPartition);

        if(readIdMap == null)
            return;

        Fragment fragment = readIdMap.get(resolvedFragment.id());

        if(fragment == null)
            return;

        if(fragment.status() == UNCLEAR)
        {
            // remove this position info as well
            Map<Integer,List<Fragment>> positionFragmentsMap = mIncompleteFragmentPositions.get(chrPartition);

            if(positionFragmentsMap != null)
            {
                // only this fragment, since the others may still turn out to be duplicates
                List<Fragment> fragments = positionFragmentsMap.get(fragment.initialPosition());

                if(fragments != null)
                {
                    fragments.remove(fragment);
                    addImpactedPosition(impactedPositions, chrPartition, fragment.initialPosition());
                }
            }
        }

        fragment.setStatus(resolvedFragment.status());
        mRecordWriter.writeFragment(fragment);
        readIdMap.remove(resolvedFragment.id());

        // could remove read ID map if empty
    }

    private void storeResolvedStatus(final String chrPartition, final Fragment resolvedFragment)
    {
        Map<String,FragmentStatus> readIdMap = mFragmentStatus.get(chrPartition);

        if(readIdMap == null)
        {
            readIdMap = Maps.newHashMap();
            mFragmentStatus.put(chrPartition, readIdMap);
        }

        readIdMap.put(resolvedFragment.id(), resolvedFragment.status());
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
            Map<String, FragmentStatus> statusMap = mFragmentStatus.get(remotePartition);

            if(statusMap != null)
            {
                FragmentStatus status = statusMap.get(supplementary.id());

                if(status != null)
                {
                    supplementary.setStatus(status);
                    mRecordWriter.writeFragment(supplementary);
                    return;
                }
            }

            // otherwise check for an unclear fragment to add this to
            Map<String,Fragment> fragmentMap = mIncompleteFragments.get(remotePartition);

            if(fragmentMap != null)
            {
                Fragment fragment = fragmentMap.get(supplementary.id());

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
            Map<String,Fragment> readIdMap = mIncompleteFragments.get(chrPartition);

            if(readIdMap == null)
            {
                readIdMap = Maps.newHashMap();
                mIncompleteFragments.put(chrPartition, readIdMap);
            }

            Fragment fragment = readIdMap.get(supplementary.id());

            if(fragment != null)
                supplementary.reads().forEach(x -> fragment.addRead(x));
            else
                readIdMap.put(supplementary.id(), supplementary);
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
                    Map<String, FragmentStatus> readIdMap = mFragmentStatus.get(remotePartition);

                    if(readIdMap == null)
                        continue;

                    FragmentStatus status = readIdMap.get(unclearFragment.id());

                    if(status != null)
                    {
                        foundMatches = true;

                        resolvedFragments.add(unclearFragment);
                        unclearFragment.setStatus(status);

                        // readIdMap.remove(unclearFragment.id()); // may be needed for other reads (eg supplementaries)
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
                Map<String,Fragment> readIdMap = mIncompleteFragments.get(remotePartition);

                if(readIdMap == null)
                    continue;

                Fragment matchedFragment = readIdMap.get(unclearFragment.id());

                if(matchedFragment != null)
                {
                    if(matchedFragment.status() == UNCLEAR)
                        foundMatches = true;

                    // add the cached fragment's reads
                    matchedFragment.reads().forEach(x -> unclearFragment.addRead(x));

                    // remove it from the remote partition and note it needs re-evaluation
                    readIdMap.remove(unclearFragment.id());

                    if(!remotePartition.equals(chrPartition))
                    {
                        // need to remove from its remote position group as well, and then deal with that
                        Map<Integer, List<Fragment>> positionFragmentsMap = mIncompleteFragmentPositions.get(remotePartition);

                        if(positionFragmentsMap != null)
                        {
                            List<Fragment> posFragments = positionFragmentsMap.get(matchedFragment.initialPosition());

                            if(posFragments != null)
                            {
                                posFragments.remove(matchedFragment);

                                // this depletes the existing position's fragment but leaves them with a status of UNCLEAR, possibly NONE
                                addImpactedPosition(impactedPositions, remotePartition, matchedFragment.initialPosition());
                            }
                        }
                    }
                    else
                    {
                        // a single-partition fragment match - another unclear or supplementaries
                    }
                }
            }
        }

        if(foundMatches)
        {
            // should be able to resolve some or all of the fragments at this position
            resolveUnclearFragments(positionFragments.Fragments);
        }

        // only store this position's fragment if there are further remote reads expected
        if(!positionFragments.Fragments.isEmpty())
        {
            Map<String,Fragment> readIdMap = mIncompleteFragments.get(chrPartition);

            if(readIdMap == null)
            {
                readIdMap = Maps.newHashMap();
                mIncompleteFragments.put(chrPartition, readIdMap);
            }

            for(Fragment fragment : positionFragments.Fragments)
            {
                readIdMap.put(fragment.id(), fragment);
            }

            Map<Integer,List<Fragment>> positionFragmentsMap = mIncompleteFragmentPositions.get(chrPartition);

            if(positionFragmentsMap == null)
            {
                positionFragmentsMap = Maps.newHashMap();
                mIncompleteFragmentPositions.put(chrPartition, positionFragmentsMap);
            }

            positionFragmentsMap.put(positionFragments.Position, positionFragments.Fragments);
        }
    }

    private void processImpactedPositions(final Map<String,Set<Integer>> impactedPositions)
    {
        for(Map.Entry<String,Set<Integer>> entry : impactedPositions.entrySet())
        {
            Map<Integer,List<Fragment>> positionFragmentsMap = mIncompleteFragmentPositions.get(entry.getKey());

            if(positionFragmentsMap == null)
                continue;

            for(Integer position : entry.getValue())
            {
                List<Fragment> positionFragments = positionFragmentsMap.get(position);

                if(positionFragments == null)
                    continue;

                resolveUnclearFragments(positionFragments);

                if(positionFragments.isEmpty())
                    positionFragmentsMap.remove(position);
            }
        }
    }

    private void resolveUnclearFragments(final List<Fragment> positionFragments)
    {
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<PositionFragments> incompletePositionFragments = Lists.newArrayList();

        classifyFragments(positionFragments, resolvedFragments, incompletePositionFragments);

        for(Fragment fragment : resolvedFragments)
        {
            positionFragments.remove(fragment);
            mRecordWriter.writeFragment(fragment);
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
            int cachedFragments = mIncompleteFragments.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.readCount()).sum()).sum();
            int cachedResolved = mFragmentStatus.values().stream().mapToInt(x -> x.size()).sum();

            BM_LOGGER.info("group cache partitions processed({}) cached fragments({}) resolved({})",
                    mProcessedPartitions.size(), cachedFragments, cachedResolved);
        }
    }

    public void handleRemaining()
    {
        int cachedFragments = mIncompleteFragments.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.readCount()).sum()).sum();
        int cachedResolved = mFragmentStatus.values().stream().mapToInt(x -> x.size()).sum();

        BM_LOGGER.info("final group cache cached fragments({}) resolved({})",
                mProcessedPartitions.size(), cachedFragments, cachedResolved);

        mIncompleteFragmentPositions.clear();

        for(Map<String,Fragment> readIdMap : mIncompleteFragments.values())
        {
            for(Fragment fragment : readIdMap.values())
            {
                fragment.setStatus(NONE);
                mRecordWriter.writeFragment(fragment);
            }
        }
    }

    public void reset()
    {
        mFragmentStatus.clear();
        mIncompleteFragmentPositions.clear();
        mIncompleteFragments.clear();
        mProcessedPartitions.clear();
    }
}
