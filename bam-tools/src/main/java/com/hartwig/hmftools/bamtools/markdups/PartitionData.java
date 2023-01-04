package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.bamtools.markdups.ResolvedFragmentState.fragmentState;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

public class PartitionData
{
    private final String mChrPartition;

    // fragment status from resolved fragments, keyed by readId
    private final Map<String,ResolvedFragmentState> mFragmentStatus;

    // supplmentary and candidate duplicate reads, keyed by readId
    private final Map<String,Fragment> mIncompleteFragments;

    // positions with candidate duplicate fragments, keyed by a unique position-based key for the group
    private final Map<String,CandidateDuplicates> mCandidateDuplicatesMap;

    // any update to the maps is done under a lock
    private Lock mLock;
    private long mLastCacheCount;
    private String mCurrentCaller;

    private static final int LOG_CACHE_COUNT = 1000;

    public PartitionData(final String chrPartition)
    {
        mChrPartition = chrPartition;
        mFragmentStatus = Maps.newHashMap();
        mIncompleteFragments = Maps.newHashMap();
        mCandidateDuplicatesMap = Maps.newHashMap();
        mLock = new ReentrantLock();
        mCurrentCaller = null;
    }

    public void processPrimaryFragments(
            final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList, final String caller)
    {
        // gather any cached mate reads, attempt to resolve any candidate duplicates and feed back the resultant set of resolved fragments
        try
        {
            mLock.lock();
            mCurrentCaller = caller;

            for(Fragment fragment : resolvedFragments)
            {
                processResolvedFragment(fragment);
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

    private void processResolvedFragment(final Fragment fragment)
    {
        if(fragment.allReadsPresent())
            return;

        // gather any higher mate or supplementary reads into this resolved fragment to be written
        Fragment existingFragment = mIncompleteFragments.get(fragment.id());

        if(existingFragment != null)
        {
            // existingFragment.setStatus(fragment.status());
            existingFragment.reads().forEach(x -> fragment.addRead(x));

            mIncompleteFragments.remove(fragment.id());

            if(fragment.allReadsPresent()) // no need to store state for reads to come
                return;
        }

        ResolvedFragmentState resolvedState = fragmentState(fragment);

        if(!resolvedState.isValid())
        {
            BM_LOGGER.error("fragment({}) invalid state({})", fragment, resolvedState);
        }

        mFragmentStatus.put(fragment.id(), resolvedState);
    }

    private boolean acquireLock()
    {
        try
        {
            boolean acquired = mLock.tryLock(200, TimeUnit.MILLISECONDS);

            if(!acquired)
            {
                BM_LOGGER.warn("partition({}) lock unacquired within time, caller({})", mChrPartition, mCurrentCaller);
            }

            return acquired;
        }
        catch(InterruptedException e)
        {
            BM_LOGGER.error(" partition({}) lock error: {}", mChrPartition, e.toString());
            e.printStackTrace();
            return false;
        }
    }

    public List<Fragment> processIncompleteFragment(final Fragment fragment, final String caller)
    {
        try
        {
            if(!acquireLock())
                return null;

            mCurrentCaller = caller;

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

                return null;
            }

            // next check for a candidate duplicate fragment to add this to
            Fragment existingFragment = mIncompleteFragments.get(fragment.id());

            if(existingFragment != null)
            {
                fragment.reads().forEach(x -> existingFragment.addRead(x));

                if(existingFragment.status() == CANDIDATE)
                {
                    // check if the set of candidates is now complete and ready for classification
                    CandidateDuplicates candidateDuplicates = mCandidateDuplicatesMap.get(existingFragment.candidateDupKey());

                    if(candidateDuplicates != null)
                    {
                        checkResolveCandidateDuplicates(candidateDuplicates);

                        if(candidateDuplicates.finalised())
                            return candidateDuplicates.fragments();
                    }
                }

                return null;
            }

            // store the new fragment
            mIncompleteFragments.put(fragment.id(), fragment);
            return null;
        }
        finally
        {
            mLock.unlock();
        }
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

    private void checkResolveCandidateDuplicates(final CandidateDuplicates candidateDuplicates)
    {
        if(!candidateDuplicates.allFragmentsReady())
            return;

        candidateDuplicates.finaliseFragmentStatus();

        for(Fragment fragment : candidateDuplicates.fragments())
        {
            mIncompleteFragments.remove(fragment.id());

            // store this new resolved state if more reads are expected for the fragment
            if(!fragment.allReadsPresent())
            {
                ResolvedFragmentState resolvedState = fragmentState(fragment);

                if(!resolvedState.isValid())
                {
                    BM_LOGGER.error("fragment({}) invalid state({})", fragment, resolvedState);
                }

                mFragmentStatus.put(fragment.id(), resolvedState);
            }
        }

        mCandidateDuplicatesMap.remove(candidateDuplicates.key());
    }

    public List<Fragment> extractRemainingFragments()
    {
        if(mIncompleteFragments.isEmpty())
            return Collections.EMPTY_LIST;

        // not under lock since called only when all partitions are complete
        BM_LOGGER.debug("final partition({}) state: incomplete({}) candidateGroups({}) resolved({})",
                mChrPartition, mIncompleteFragments.size(), mCandidateDuplicatesMap.size(), mFragmentStatus.size());

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

        BM_LOGGER.debug("partition({}) state: incomplete({}) candidateGroups({}) resolved({})",
                mChrPartition, mIncompleteFragments.size(), mCandidateDuplicatesMap.size(), mFragmentStatus.size());
    }

    public void logCacheCounts()
    {
        try
        {
            mLock.lock();

            BM_LOGGER.debug("partition({}) state: incomplete({}) candidateGroups({}) resolved({})",
                    mChrPartition, mIncompleteFragments.size(), mCandidateDuplicatesMap.size(), mFragmentStatus.size());
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
    public List<Fragment> processIncompleteFragment(final Fragment fragment)
    {
        return processIncompleteFragment(fragment, null);
    }

    @VisibleForTesting
    public void processPrimaryFragments(final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList)
    {
        processPrimaryFragments(resolvedFragments, candidateDuplicatesList, null);
    }
}
