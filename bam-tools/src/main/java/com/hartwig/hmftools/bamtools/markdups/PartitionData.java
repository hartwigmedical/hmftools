package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.ResolvedFragmentState.fragmentState;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

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

    private Lock mLock;

    private static int MAX_FRAGMENT_READS = 4;

    public PartitionData(final String chrPartition)
    {
        mChrPartition = chrPartition;
        mFragmentStatus = Maps.newHashMap();
        mIncompleteFragments = Maps.newHashMap();
        mCandidateDuplicatesMap = Maps.newHashMap();
        mLock = new ReentrantLock();
    }

    public int incompleteFragments()
    {
        return mIncompleteFragments.size();
    }
    public int resolvedFragments()
    {
        return mFragmentStatus.size();
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

        int processedReads = fragment.readCount();
        int expectedReads = fragment.expectedReadCount();

        if(processedReads >= expectedReads)
        {
            BM_LOGGER.error("fragment({}) invalid state: reads expected({}) processed({})",
                    fragment, expectedReads, processedReads);
            return;
        }

        ResolvedFragmentState resolvedState = fragmentState(fragment);
        mFragmentStatus.put(fragment.id(), resolvedState);
    }

    public List<Fragment> processIncompleteFragment(final Fragment fragment)
    {
        // a supplementary or higher mate read - returns any resolved fragments resulting from add this new read
        List<Fragment> resolvedFragments = null;

        try
        {
            mLock.lock();

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
                ResolvedFragmentState resolvedState = fragmentState(fragment);;
                mFragmentStatus.put(fragment.id(), resolvedState);
            }
        }

        return resolvedFragments;
    }

    /*
    public List<Fragment> extractResolvedFragments()
    {
        try
        {
            List<Fragment> resolvedFragments = Lists.newArrayList(mResolvedFragments);
            mResolvedFragments.clear();

            resolvedFragments.addAll(findResolvedFragments());
            return resolvedFragments;
        }
        finally
        {
            mLock.unlock();
        }
    }
    */

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

    public String toString()
    {
        return format("%s status(%d) frags(%d) positions(%d)",
            mChrPartition, mFragmentStatus.size(), mIncompleteFragments.size(), mCandidateDuplicatesMap.size());
    }

    @VisibleForTesting
    public Map<String,ResolvedFragmentState> fragmentStatusMap() { return mFragmentStatus; }

    @VisibleForTesting
    public Map<String,Fragment> incompleteFragmentMap() { return mIncompleteFragments; }

    @VisibleForTesting
    public Map<Integer,List<Fragment>> candidateDuplicatesMap() { return mCandidateDuplicatesMap; }
}
