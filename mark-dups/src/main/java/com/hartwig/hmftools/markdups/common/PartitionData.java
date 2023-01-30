package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.markdups.common.ResolvedFragmentState.fragmentState;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.markdups.umi.UmiConfig;
import com.hartwig.hmftools.markdups.umi.UmiGroup;

public class PartitionData
{
    private final String mChrPartition;

    // fragment status from resolved fragments, keyed by readId
    private final Map<String,ResolvedFragmentState> mFragmentStatus;

    private final Map<String,UmiGroup> mUmiGroups;

    // supplmentary and candidate duplicate reads, keyed by readId
    private final Map<String,Fragment> mIncompleteFragments;

    // positions with candidate duplicate fragments, keyed by a unique position-based key for the group
    private final Map<String,CandidateDuplicates> mCandidateDuplicatesMap;

    private final DuplicateGroups mDuplicateGroups;

    // any update to the maps is done under a lock
    private Lock mLock;
    private long mLastCacheCount;

    private Set<UmiGroup> mUpdatedUmiGroups;
    private Set<CandidateDuplicates> mUpdatedCandidateDuplicates;

    private static final int LOG_CACHE_COUNT = 1000;

    public PartitionData(final String chrPartition, final UmiConfig umiConfig)
    {
        mChrPartition = chrPartition;
        mFragmentStatus = Maps.newHashMap();
        mIncompleteFragments = Maps.newHashMap();
        mCandidateDuplicatesMap = Maps.newHashMap();
        mUmiGroups = Maps.newHashMap();
        mDuplicateGroups = new DuplicateGroups(umiConfig);
        mUpdatedUmiGroups = Sets.newHashSet();
        mUpdatedCandidateDuplicates = Sets.newHashSet();

        mLock = new ReentrantLock();
    }

    public Statistics statistics() { return mDuplicateGroups.statistics(); }

    public void processPrimaryFragments(
            final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList, final List<UmiGroup> umiGroups)
    {
        // gather any cached mate reads, attempt to resolve any candidate duplicates and feed back the resultant set of resolved fragments
        try
        {
            mLock.lock();

            resolvedFragments.forEach(x -> processResolvedFragment(x));

            if(umiGroups != null)
            {
                umiGroups.forEach(x -> processUmiGroup(x));
            }

            for(CandidateDuplicates candidateDuplicates : candidateDuplicatesList)
            {
                processCandidateDuplicates(candidateDuplicates);

                // add any additional resolved fragments after gathering mate reads
                if(candidateDuplicates.finalised())
                    resolvedFragments.addAll(candidateDuplicates.fragments());
            }
        }
        finally
        {
            checkCachedCounts();
            mLock.unlock();
        }
    }

    private void processUmiGroup(final UmiGroup umiGroup)
    {
        if(umiGroup.allReadsReceived())
            return;

        boolean addedRead = false;

        // for(Fragment fragment : umiGroup.fragments())
        List<String> readIds = umiGroup.getReadIds();
        for(String readId : readIds)
        {
            Fragment existingFragment = mIncompleteFragments.get(readId);

            if(existingFragment != null)
            {
                existingFragment.reads().forEach(x -> umiGroup.addRead(x));

                mIncompleteFragments.remove(readId);

                addedRead = true;
            }
        }

        if(addedRead && umiGroup.allReadsReceived())
            return;

        // store the UMI group to pick up mates and supplementaries when they arrive
        for(String readId : readIds)
        {
            mUmiGroups.put(readId, umiGroup);
        }
    }

    private void processResolvedFragment(final Fragment fragment)
    {
        if(fragment.allReadsPresent())
            return;

        // gather any higher mate or supplementary reads into this resolved fragment to be written
        Fragment existingFragment = mIncompleteFragments.get(fragment.id());

        if(existingFragment != null)
        {
            existingFragment.reads().forEach(x -> fragment.addRead(x));

            mIncompleteFragments.remove(fragment.id());

            if(fragment.allReadsPresent()) // no need to store state for reads to come
                return;
        }

        if(fragment.status() != NONE && umiEnabled())
            return;

        ResolvedFragmentState resolvedState = fragmentState(fragment);

        if(!resolvedState.isValid())
        {
            MD_LOGGER.error("fragment({}) invalid state({})", fragment, resolvedState);
        }

        mFragmentStatus.put(fragment.id(), resolvedState);
    }

    private void processCandidateDuplicates(final CandidateDuplicates candidateDuplicates)
    {
        // this position cannot already exist
        // the resolved status cannot be known since this contains the lower reads
        // there may be mate reads and supplementaries to add to the fragments it contains
        boolean hasCompleteReads = false;

        for(Fragment fragment : candidateDuplicates.fragments())
        {
            Fragment existingFragment = mIncompleteFragments.get(fragment.id());

            if(existingFragment != null)
            {
                existingFragment.reads().forEach(x -> fragment.addRead(x));
                mIncompleteFragments.put(fragment.id(), fragment); // replace it

                if(existingFragment.primaryReadsPresent())
                    hasCompleteReads = true;
            }
        }

        if(hasCompleteReads)
        {
            checkResolveCandidateDuplicates(candidateDuplicates);
        }

        if(!candidateDuplicates.finalised())
        {
            for(Fragment fragment : candidateDuplicates.fragments())
            {
                if(!mIncompleteFragments.containsValue(fragment.id()))
                    mIncompleteFragments.put(fragment.id(), fragment);
            }

            mCandidateDuplicatesMap.put(candidateDuplicates.key(), candidateDuplicates);
        }
    }

    public PartitionResults processIncompleteFragments(final List<Fragment> fragments)
    {
        try
        {
            mLock.lock();

            PartitionResults partitionResults = new PartitionResults();

            for(Fragment fragment : fragments)
            {
                handleIncompleteFragment(fragment);

                if(fragment.status().isResolved())
                    partitionResults.addResolvedFragment(fragment);
            }

            processUpdatedGroups(partitionResults);

            return partitionResults;
        }
        finally
        {
            mLock.unlock();
        }
    }

    public PartitionResults processIncompleteFragment(final Fragment fragment)
    {
        try
        {
            mLock.lock();
            boolean matched = handleIncompleteFragment(fragment);

            if(!matched)
                return null;

            // only create results if the fragment is part of a group
            if(fragment.status() == NONE)
                return null;

            PartitionResults partitionResults = new PartitionResults();

            processUpdatedGroups(partitionResults);
            return partitionResults;
        }
        finally
        {
            mLock.unlock();
        }
    }

    private boolean handleIncompleteFragment(final Fragment fragment)
    {
        // a supplementary or higher mate read - returns any resolved fragments resulting from add this new read

        // first look for a resolved status
        ResolvedFragmentState resolvedState = mFragmentStatus.get(fragment.id());

        if(resolvedState != null)
        {
            fragment.setStatus(resolvedState.Status);

            // update the resolved state
            resolvedState.update(fragment.reads());

            if(resolvedState.allReceived())
                mFragmentStatus.remove(fragment.id());

            return true;
        }

        UmiGroup umiGroup = mUmiGroups.get(fragment.id());

        if(umiGroup != null)
        {
            fragment.reads().forEach(x -> umiGroup.addRead(x));
            mUpdatedUmiGroups.add(umiGroup);
            return true;
        }

        // next check for a UMI group or candidate duplicate group to add this to
        Fragment existingFragment = mIncompleteFragments.get(fragment.id());

        if(existingFragment != null)
        {
            fragment.reads().forEach(x -> existingFragment.addRead(x));

            if(existingFragment.status() == CANDIDATE && existingFragment.primaryReadsPresent())
            {
                // check if the set of candidates is now complete and ready for classification
                CandidateDuplicates candidateDuplicates = mCandidateDuplicatesMap.get(existingFragment.candidateDupKey());

                if(candidateDuplicates != null)
                {
                    mUpdatedCandidateDuplicates.add(candidateDuplicates);
                    return true;
                }
            }

            return false;
        }

        // store the new fragment
        mIncompleteFragments.put(fragment.id(), fragment);
        return false;
    }

    private void storeUmiGroup(final UmiGroup umiGroup)
    {
        if(umiGroup.allReadsReceived())
            return;

        umiGroup.fragments().forEach(x -> mUmiGroups.put(x.id(), umiGroup));
    }

    private void checkRemoveUmiGroup(final UmiGroup umiGroup)
    {
        if(!umiGroup.allReadsReceived())
            return;

        // remove by each read ID
        List<String> groupReadIds = umiGroup.getReadIds();

        if(groupReadIds == null)
        {
            MD_LOGGER.error("umiGroup({}) has no read IDs: {}", umiGroup.umiId(), umiGroup.toString());
            return;
        }

        groupReadIds.forEach(x -> mUmiGroups.remove(x));
    }

    private void checkResolveCandidateDuplicates(final CandidateDuplicates candidateDuplicates)
    {
        if(!candidateDuplicates.allFragmentsReady())
            return;

        List<List<Fragment>> duplicateGroups = candidateDuplicates.finaliseFragmentStatus();

        if(umiEnabled())
        {
            List<UmiGroup> umiGroups = mDuplicateGroups.processDuplicateUmiGroups(duplicateGroups);

            for(UmiGroup umiGroup : umiGroups)
            {
                mUpdatedUmiGroups.add(umiGroup);

                // store only if incomplete
                storeUmiGroup(umiGroup);
            }
        }
        else
        {
            mDuplicateGroups.processDuplicateGroups(duplicateGroups);
        }

        for(Fragment fragment : candidateDuplicates.fragments())
        {
            mIncompleteFragments.remove(fragment.id());

            // store this new resolved state if more reads are expected for the fragment
            if(fragment.allReadsPresent())
                continue;

            if(fragment.umiId() != null) // cached with the UMI group
                continue;

            ResolvedFragmentState resolvedState = fragmentState(fragment);

            if(!resolvedState.isValid())
            {
                MD_LOGGER.error("fragment({}) invalid state({})", fragment, resolvedState);
            }

            mFragmentStatus.put(fragment.id(), resolvedState);
        }

        mCandidateDuplicatesMap.remove(candidateDuplicates.key());
    }

    private void processUpdatedGroups(final PartitionResults partitionResults)
    {
        if(mUpdatedUmiGroups.isEmpty() && mUpdatedCandidateDuplicates.isEmpty())
            return;

        for(CandidateDuplicates candidateDuplicates : mUpdatedCandidateDuplicates)
        {
            checkResolveCandidateDuplicates(candidateDuplicates);

            if(candidateDuplicates.finalised())
                partitionResults.addResolvedFragments(candidateDuplicates.fragments());
        }

        mUpdatedCandidateDuplicates.clear();

        // only add UMI groups if they have complete sets of reads

        for(UmiGroup umiGroup : mUpdatedUmiGroups)
        {
            if(umiGroup.hasCompleteReadGroup())
            {
                partitionResults.addUmiGroup(umiGroup);
                checkRemoveUmiGroup(umiGroup);
            }
        }

        mUpdatedUmiGroups.clear();
    }

    private boolean umiEnabled() { return mDuplicateGroups.umiConfig().Enabled; }

    public List<Fragment> extractRemainingFragments()
    {
        if(mIncompleteFragments.isEmpty())
            return Collections.EMPTY_LIST;

        // not under lock since called only when all partitions are complete
        MD_LOGGER.debug("partition({}) final state: {}", mChrPartition, cacheCountsStr());

        List<Fragment> remainingFragments = mIncompleteFragments.values().stream().collect(Collectors.toList());

        mFragmentStatus.clear();
        mIncompleteFragments.clear();
        mCandidateDuplicatesMap.clear();

        return remainingFragments;
    }

    private void checkCachedCounts()
    {
        long cacheCount = mIncompleteFragments.size() + mFragmentStatus.size();

        if(abs(mLastCacheCount - cacheCount) < LOG_CACHE_COUNT)
            return;

        mLastCacheCount = cacheCount;

        MD_LOGGER.debug("partition({}) check state: {}}", mChrPartition, cacheCountsStr());
    }

    private String cacheCountsStr()
    {
        long incompleteSupp = mIncompleteFragments.values().stream().filter(x -> x.status() == SUPPLEMENTARY).count();
        int maxCandidateGroup = mCandidateDuplicatesMap.values().stream().mapToInt(x -> x.fragmentCount()).max().orElse(0);
        long resolvedNoSupp = mFragmentStatus.values().stream().filter(x -> x.MateReceived).count();
        long resolvedNoMate = mFragmentStatus.values().stream().filter(x -> x.ProcessedSupplementaries < x.ExpectedSupplementaries).count();

        return format("incomplete(%d supp=%d) candidateGroups(%d max=%d) resolved(%d supp=%d mate=%d) umiGroupReads(%d)",
                mIncompleteFragments.size(), incompleteSupp, mCandidateDuplicatesMap.size(), maxCandidateGroup,
                mFragmentStatus.size(), resolvedNoSupp, resolvedNoMate, mUmiGroups.size());
    }

    public void logCacheCounts()
    {
        try
        {
            mLock.lock();

            MD_LOGGER.debug("partition({}) log state: {}", mChrPartition, cacheCountsStr());
        }
        finally
        {
            mLock.unlock();
        }
    }

    public void clear()
    {
        try
        {
            mLock.lock();

            mFragmentStatus.clear();
            mIncompleteFragments.clear();
            mCandidateDuplicatesMap.clear();
            mUmiGroups.clear();
        }
        finally
        {
            mLock.unlock();
        }
    }

    @VisibleForTesting
    public Map<String,ResolvedFragmentState> fragmentStatusMap() { return mFragmentStatus; }

    @VisibleForTesting
    public Map<String,Fragment> incompleteFragmentMap() { return mIncompleteFragments; }

    @VisibleForTesting
    public Map<String,CandidateDuplicates> candidateDuplicatesMap() { return mCandidateDuplicatesMap; }

    @VisibleForTesting
    public Map<String,ResolvedFragmentState> resolvedFragmentStateMap() { return mFragmentStatus; }

    @VisibleForTesting
    public Map<String,UmiGroup> umiGroupMap() { return mUmiGroups; }

    @VisibleForTesting
    public void processPrimaryFragments(final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList)
    {
        processPrimaryFragments(resolvedFragments, candidateDuplicatesList, Collections.EMPTY_LIST);
    }
}
