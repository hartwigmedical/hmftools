package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.ResolvedFragmentState.fragmentState;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

public class PartitionData
{
    private final String mChrPartition;

    // fragment status from resolved fragments, keyed by readId
    private final Map<String,ResolvedFragmentState> mFragmentStatus;

    // supplmentary and candidate duplicate reads, keyed by readId
    private final Map<String,Fragment> mIncompleteFragments;

    // positions with candidate duplicate fragments, keyed by initial fragment coordinate position
    private final Map<Integer,List<Fragment>> mCandidateDuplicatesMap;

    // any update to the maps is done under a lock
    private Lock mLock;

    public PartitionData(final String chrPartition)
    {
        mChrPartition = chrPartition;
        mFragmentStatus = Maps.newHashMap();
        mIncompleteFragments = Maps.newHashMap();
        mCandidateDuplicatesMap = Maps.newHashMap();
        mLock = new ReentrantLock();
    }

    public void processPrimaryFragments(final List<Fragment> resolvedFragments, @Nullable final CandidateDuplicates candidateDuplicates)
    {
        // gather any cached mate reads, attempt to resolve any candidate duplicates and feed back the resultant set of resolved fragments
        try
        {
            mLock.lock();

            for(Fragment fragment : resolvedFragments)
            {
                processResolvedFragment(fragment);
            }

            if(candidateDuplicates != null)
            {
                List<Fragment> resolvedCandidates = processCandidateDuplicates(candidateDuplicates);

                // add any additional resolved fragments after gathering mate reads
                if(resolvedCandidates != null)
                    resolvedFragments.addAll(resolvedCandidates);
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
            return mLock.tryLock(2, TimeUnit.SECONDS);
        }
        catch(InterruptedException e)
        {
            BM_LOGGER.error(" partition({}) lock error: {}", mChrPartition, e.toString());
            e.printStackTrace();
            return false;
        }
    }

    public List<Fragment> processIncompleteFragment(final Fragment fragment)
    {
        try
        {
            if(!acquireLock())
            {
                BM_LOGGER.warn("partition({}) process primary failed to acquire lock within time", mChrPartition);
                return null;
            }

            // a supplementary or higher mate read - returns any resolved fragments resulting from add this new read
            List<Fragment> resolvedFragments = null;

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

                if(existingFragment.status() == UNCLEAR)
                {
                    // check for a now complete group
                    resolvedFragments = findResolvedFragments(existingFragment);
                }

                return resolvedFragments;
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

    private List<Fragment> findResolvedFragments(final Fragment candidateDuplicate)
    {
        List<Fragment> candidateFragments = mCandidateDuplicatesMap.get(candidateDuplicate.initialPosition());

        if(candidateFragments == null)
            return null;

        List<Fragment> resolvedFragments = resolveCandidateDuplicates(candidateFragments);

        if(candidateFragments.isEmpty())
            mCandidateDuplicatesMap.remove(candidateDuplicate.initialPosition());

        return resolvedFragments;
    }

    private List<Fragment> processCandidateDuplicates(final CandidateDuplicates candidateDuplicates)
    {
        // this position cannot already exist
        // the resolved status cannot be known since this contains the lower reads
        // there may be mate reads and supplementaries to add to the fragments it contains
        boolean hasCompleteReads = false;

        for(Fragment fragment : candidateDuplicates.Fragments)
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

        List<Fragment> resolvedFragments = null;

        if(hasCompleteReads)
        {
            resolvedFragments = resolveCandidateDuplicates(candidateDuplicates.Fragments);
        }

        if(!candidateDuplicates.Fragments.isEmpty())
        {
            for(Fragment fragment : candidateDuplicates.Fragments)
            {
                if(!mIncompleteFragments.containsValue(fragment.id()))
                    mIncompleteFragments.put(fragment.id(), fragment);
            }

            mCandidateDuplicatesMap.put(candidateDuplicates.Position, candidateDuplicates.Fragments);
        }

        return resolvedFragments;
    }

    private List<Fragment> resolveCandidateDuplicates(final List<Fragment> positionFragments)
    {
        if(positionFragments.isEmpty())
            return null;

        List<Fragment> resolvedFragments = Lists.newArrayList();

        positionFragments.forEach(x -> x.setStatus(UNSET)); // reset before evaluation

        classifyFragments(positionFragments, resolvedFragments, null);

        for(Fragment fragment : resolvedFragments)
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

        return resolvedFragments;
    }

    public List<Fragment> extractRemainingFragments()
    {
        if(mIncompleteFragments.isEmpty())
            return Collections.EMPTY_LIST;

        // not under lock since called only when all partitions are complete
        BM_LOGGER.debug("final partition({}) state: incomplete({}) positions({}) resolved({})",
                mChrPartition, mIncompleteFragments.size(), mCandidateDuplicatesMap.size(), mFragmentStatus.size());

        List<Fragment> remainingFragments = mIncompleteFragments.values().stream().collect(Collectors.toList());

        mFragmentStatus.clear();
        mIncompleteFragments.clear();
        mCandidateDuplicatesMap.clear();

        return remainingFragments;
    }

    private static final int LOG_CACHE_COUNT = 10000;
    private long mLastCacheCount;

    private void checkCachedCounts()
    {
        if(mLastCacheCount < LOG_CACHE_COUNT)
            return;

        long cacheCount = mIncompleteFragments.size() + mFragmentStatus.size();

        if(abs(mLastCacheCount - cacheCount) < LOG_CACHE_COUNT)
            return;

        mLastCacheCount = cacheCount;

        BM_LOGGER.debug("partition({}) state: incomplete({}) positions({}) resolved({})",
                mChrPartition, mIncompleteFragments.size(), mCandidateDuplicatesMap.size(), mFragmentStatus.size());
    }

    public void logCacheCounts()
    {
        try
        {
            mLock.lock();

            BM_LOGGER.debug("partition({}) state: incomplete({}) positions({}) resolved({})",
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
    public Map<Integer,List<Fragment>> candidateDuplicatesMap() { return mCandidateDuplicatesMap; }

        /*
    public String toString()
    {
        return format("%s status(%d) frags(%d) positions(%d)",
            mChrPartition, mFragmentStatus.size(), mIncompleteFragments.size(), mCandidateDuplicatesMap.size());
    }
    */
}
