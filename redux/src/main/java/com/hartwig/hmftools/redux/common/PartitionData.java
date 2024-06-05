package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.NANO_IN_MILLISECOND;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.redux.common.FragmentUtils.readToString;
import static com.hartwig.hmftools.redux.common.ReadMatch.NO_READ_MATCH;
import static com.hartwig.hmftools.redux.common.ResolvedFragmentState.fragmentState;

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
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.write.BamWriter;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public class PartitionData
{
    private final String mChrPartition;

    // fragment status from resolved fragments, keyed by readId
    private final Map<String,ResolvedFragmentState> mFragmentStatus;

    private final Map<String,DuplicateGroup> mDuplicateGroupMap; // keyed by readId

    // supplmentary and candidate duplicate reads, keyed by readId
    private final Map<String,Fragment> mIncompleteFragments;

    // positions with candidate duplicate fragments, keyed by a unique position-based key for the group
    private final Map<String,CandidateDuplicates> mCandidateDuplicatesMap;

    private final DuplicateGroupBuilder mDuplicateGroupBuilder;

    // any update to the maps is done under a lock
    private Lock mLock;
    private long mLastCacheCount;
    private long mLockAcquireTime;
    private boolean mPerfChecks;

    private Set<DuplicateGroup> mUpdatedDuplicateGroups;
    private Set<CandidateDuplicates> mUpdatedCandidateDuplicates;

    private static final int LOG_CACHE_COUNT = 50000;

    public PartitionData(final String chrPartition, final ReduxConfig config)
    {
        mChrPartition = chrPartition;
        mFragmentStatus = Maps.newHashMap();
        mIncompleteFragments = Maps.newHashMap();
        mCandidateDuplicatesMap = Maps.newHashMap();
        mDuplicateGroupMap = Maps.newHashMap();
        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mUpdatedDuplicateGroups = Sets.newHashSet();
        mUpdatedCandidateDuplicates = Sets.newHashSet();

        mLock = new ReentrantLock();
        mLockAcquireTime = 0;
        mPerfChecks = false;
    }

    public String partitionStr() { return mChrPartition; }
    public Statistics statistics() { return mDuplicateGroupBuilder.statistics(); }

    public void togglePerfChecks() { mPerfChecks = true; }
    public double totalLockTimeMs() { return mLockAcquireTime / NANO_IN_MILLISECOND; }

    public void processPrimaryFragments(
            final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList, final List<DuplicateGroup> duplicateGroups)
    {
        // gather any cached mate reads, attempt to resolve any candidate duplicates and feed back the resultant set of resolved fragments
        try
        {
            acquireLock();

            // UMIs are filtered here since the fragments themselves don't need to collect incomplete reads nor set resolved status
            resolvedFragments.stream().filter(x -> x.umi() == null).forEach(x -> processResolvedFragment(x));

            if(duplicateGroups != null)
            {
                duplicateGroups.forEach(x -> processDuplicateGroup(x));
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

    private void processDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        if(duplicateGroup.allReadsReceived())
            return;

        boolean addedRead = false;

        List<String> readIds = duplicateGroup.getReadIds();

        for(String readId : readIds)
        {
            Fragment existingFragment = mIncompleteFragments.get(readId);

            if(existingFragment != null)
            {
                existingFragment.reads().forEach(x -> duplicateGroup.addRead(x));

                mIncompleteFragments.remove(readId);
                addedRead = true;
            }
        }

        if(addedRead && duplicateGroup.allReadsReceived())
            return;

        // store the UMI group to pick up mates and supplementaries when they arrive
        for(String readId : readIds)
        {
            mDuplicateGroupMap.put(readId, duplicateGroup);
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

        // if(fragment.status() != NONE && umiEnabled()) // fragments in UMI groups don't get through this routine
        //    return;

        ResolvedFragmentState resolvedState = fragmentState(fragment);
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

    public PartitionResults processIncompleteFragments(final List<SAMRecord> reads)
    {
        try
        {
            acquireLock();

            PartitionResults partitionResults = new PartitionResults();

            for(SAMRecord read : reads)
            {
                ReadMatch readMatch = handleIncompleteFragment(read);

                if(readMatch.Status != null && readMatch.Status.isResolved())
                {
                    Fragment fragment = new Fragment(read);
                    fragment.setStatus(readMatch.Status);
                    partitionResults.addResolvedFragment(fragment);
                }
            }

            processUpdatedGroups(partitionResults);

            return partitionResults;
        }
        finally
        {
            mLock.unlock();
        }
    }

    public PartitionResults processIncompleteFragment(final SAMRecord read)
    {
        try
        {
            acquireLock();
            ReadMatch readMatch = handleIncompleteFragment(read);

            if(!readMatch.Matched)
                return null;

            // only create results if the fragment is part of a group
            PartitionResults partitionResults = new PartitionResults();

            if(readMatch.Status != null)
                partitionResults.setFragmentStatus(readMatch.Status);

            if(readMatch.Status == null || readMatch.Status != NONE)
                processUpdatedGroups(partitionResults);

            return partitionResults;
        }
        finally
        {
            mLock.unlock();
        }
    }

    private ReadMatch handleIncompleteFragment(final SAMRecord read)
    {
        // a supplementary or higher mate read - returns any resolved fragments resulting from add this new read

        // first look for a resolved status
        ResolvedFragmentState resolvedState = mFragmentStatus.get(read.getReadName());

        if(resolvedState != null)
        {
            resolvedState.update(read);

            if(resolvedState.allReceived())
                mFragmentStatus.remove(read.getReadName());

            return new ReadMatch(true, resolvedState.Status);
        }

        DuplicateGroup duplicateGroup = mDuplicateGroupMap.get(read.getReadName());

        if(duplicateGroup != null)
        {
            duplicateGroup.addRead(read);
            mUpdatedDuplicateGroups.add(duplicateGroup);
            return new ReadMatch(true, null);
        }

        // next check for a UMI group or candidate duplicate group to add this to
        Fragment existingFragment = mIncompleteFragments.get(read.getReadName());

        if(existingFragment != null)
        {
            existingFragment.addRead(read);

            if(existingFragment.status() == CANDIDATE && existingFragment.primaryReadsPresent())
            {
                // check if the set of candidates is now complete and ready for classification
                CandidateDuplicates candidateDuplicates = mCandidateDuplicatesMap.get(existingFragment.candidateDupKey());

                if(candidateDuplicates != null)
                {
                    mUpdatedCandidateDuplicates.add(candidateDuplicates);
                    return new ReadMatch(true, null);
                }
            }

            return NO_READ_MATCH;
        }

        // store the new fragment
        mIncompleteFragments.put(read.getReadName(), new Fragment(read));
        return NO_READ_MATCH;
    }

    private void storeDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        if(duplicateGroup.allReadsReceived())
            return;

        duplicateGroup.fragments().forEach(x -> mDuplicateGroupMap.put(x.id(), duplicateGroup));
    }

    private void checkRemoveUmiGroup(final DuplicateGroup duplicateGroup)
    {
        if(!duplicateGroup.allReadsReceived())
            return;

        // remove by each read ID
        List<String> groupReadIds = duplicateGroup.getReadIds();
        groupReadIds.forEach(x -> mDuplicateGroupMap.remove(x));
    }

    private void checkResolveCandidateDuplicates(final CandidateDuplicates candidateDuplicates)
    {
        if(!candidateDuplicates.allFragmentsReady())
            return;

        List<List<Fragment>> rawDuplicateGroups = candidateDuplicates.finaliseFragmentStatus(mDuplicateGroupBuilder.umiConfig().Enabled);

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                rawDuplicateGroups, false, Collections.EMPTY_LIST);

        if(duplicateGroups != null)
        {
            for(DuplicateGroup duplicateGroup : duplicateGroups)
            {
                mUpdatedDuplicateGroups.add(duplicateGroup);

                // store only if incomplete
                storeDuplicateGroup(duplicateGroup);
            }
        }

        for(Fragment fragment : candidateDuplicates.fragments())
        {
            mIncompleteFragments.remove(fragment.id());

            // store this new resolved state if more reads are expected for the fragment
            if(fragment.allReadsPresent())
                continue;

            if(fragment.umi() != null) // if not null, then part of a duplicate group and status is derived from the group
                continue;

            ResolvedFragmentState resolvedState = fragmentState(fragment);
            mFragmentStatus.put(fragment.id(), resolvedState);
        }

        mCandidateDuplicatesMap.remove(candidateDuplicates.key());
    }

    private void processUpdatedGroups(final PartitionResults partitionResults)
    {
        if(mUpdatedDuplicateGroups.isEmpty() && mUpdatedCandidateDuplicates.isEmpty())
            return;

        for(CandidateDuplicates candidateDuplicates : mUpdatedCandidateDuplicates)
        {
            checkResolveCandidateDuplicates(candidateDuplicates);

            if(candidateDuplicates.finalised())
                partitionResults.addResolvedFragments(candidateDuplicates.fragments());
        }

        mUpdatedCandidateDuplicates.clear();

        // only add UMI groups if they have complete sets of reads

        for(DuplicateGroup duplicateGroup : mUpdatedDuplicateGroups)
        {
            if(duplicateGroup.hasCompleteReadGroup())
            {
                partitionResults.addUmiGroup(duplicateGroup);
                checkRemoveUmiGroup(duplicateGroup);
            }
        }

        mUpdatedDuplicateGroups.clear();
    }

    public int writeRemainingReads(final BamWriter recordWriter, final ConsensusReads consensusReads, boolean logCachedReads)
    {
        if(mDuplicateGroupMap.isEmpty() && mIncompleteFragments.isEmpty() && mFragmentStatus.isEmpty())
            return 0;

        if(logCachedReads && RD_LOGGER.isDebugEnabled())
        {
            // not under lock since called only when all partitions are complete
            RD_LOGGER.debug("partition({}) final state: {}", mChrPartition, cacheCountsStr());
        }

        // a clean-up routine for cached fragments
        Set<DuplicateGroup> processedDuplicateGroups = Sets.newHashSet();

        int cachedReadCount = 0;

        for(DuplicateGroup duplicateGroup : mDuplicateGroupMap.values())
        {
            if(processedDuplicateGroups.contains(duplicateGroup))
                continue;

            processedDuplicateGroups.add(duplicateGroup);

            for(String readId : duplicateGroup.getReadIds())
            {
                Fragment incompleteFragment = mIncompleteFragments.get(readId);

                if(incompleteFragment != null)
                {
                    mIncompleteFragments.remove(readId);
                    incompleteFragment.reads().forEach(x -> duplicateGroup.addRead(x));
                }
            }

            int cachedUmiReads = duplicateGroup.cachedReadCount();

            if(cachedUmiReads == 0)
                continue;

            cachedReadCount += cachedUmiReads;

            List<SAMRecord> completeReads = duplicateGroup.popCompletedReads(consensusReads, true);
            recordWriter.writeDuplicateGroup(duplicateGroup, completeReads);

            if(logCachedReads)
            {
                RD_LOGGER.debug("writing {} cached reads for umi group({}) coords({})",
                        cachedUmiReads, duplicateGroup.toString(), duplicateGroup.coordinatesKey());

                for(SAMRecord read : completeReads)
                {
                    if(read.getSupplementaryAlignmentFlag() || read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                        continue;

                    RD_LOGGER.debug("writing umi read: {}", readToString(read));
                }
            }
        }

        for(Fragment fragment : mIncompleteFragments.values())
        {
            recordWriter.writeFragment(fragment);

            if(logCachedReads)
            {
                for(SAMRecord read : fragment.reads())
                {
                    if(read.getSupplementaryAlignmentFlag())
                        continue;

                    ++cachedReadCount;
                    RD_LOGGER.debug("writing incomplete read: {} status({})", readToString(read), fragment.status());
                }
            }
        }

        if(logCachedReads && !mFragmentStatus.isEmpty())
        {
            for(Map.Entry<String,ResolvedFragmentState> entry : mFragmentStatus.entrySet())
            {
                RD_LOGGER.debug("cached resolved status: {} : {}", entry.getKey(), entry.getValue());
            }
        }

        mFragmentStatus.clear();
        mIncompleteFragments.clear();
        mCandidateDuplicatesMap.clear();
        mDuplicateGroupMap.clear();

        return cachedReadCount;
    }

    private void checkCachedCounts()
    {
        long cacheCount = mIncompleteFragments.size() + mFragmentStatus.size();

        if(abs(mLastCacheCount - cacheCount) < LOG_CACHE_COUNT)
            return;

        mLastCacheCount = cacheCount;

        RD_LOGGER.debug("partition({}) check state: {}}", mChrPartition, cacheCountsStr());
    }

    private String cacheCountsStr()
    {
        long incompleteSupp = mIncompleteFragments.values().stream().filter(x -> x.status() == SUPPLEMENTARY).count();
        int maxCandidateGroup = mCandidateDuplicatesMap.values().stream().mapToInt(x -> x.fragmentCount()).max().orElse(0);
        long resolvedNoSupp = mFragmentStatus.values().stream().filter(x -> x.MateReceived).count();
        long resolvedNoMate = mFragmentStatus.values().stream().filter(x -> x.ProcessedSupplementaries < x.ExpectedSupplementaries).count();

        long umiReads = 0;
        Set<DuplicateGroup> uniqueGroups = Sets.newHashSet();

        if(mPerfChecks)
        {
            mDuplicateGroupMap.values().forEach(x -> uniqueGroups.add(x));
            umiReads = uniqueGroups.stream().mapToInt(x -> x.cachedReadCount()).sum();
        }

        return format("incomplete(%d supp=%d) resolved(%d supp=%d mate=%d) umi(groups=%d frags=%s reads=%d) candidateGroups(%d max=%d)",
                mIncompleteFragments.size(), incompleteSupp, mFragmentStatus.size(), resolvedNoSupp, resolvedNoMate,
                uniqueGroups.size(), mDuplicateGroupMap.size(), umiReads, mCandidateDuplicatesMap.size(), maxCandidateGroup);
    }

    public void logCacheCounts()
    {
        try
        {
            acquireLock();

            RD_LOGGER.debug("partition({}) log state: {}", mChrPartition, cacheCountsStr());
        }
        finally
        {
            mLock.unlock();
        }
    }

    private void acquireLock()
    {
        if(!mPerfChecks)
        {
            mLock.lock();
            return;
        }

        long startTime = System.nanoTime();
        mLock.lock();
        mLockAcquireTime += System.nanoTime() - startTime;
    }

    public String toString()
    {
        return format("%s: status(%d) incomplete(%d) candidates(%d) umis(%d)",
                mChrPartition, mFragmentStatus.size(), mIncompleteFragments.size(), mCandidateDuplicatesMap.size(), mDuplicateGroupMap.size());
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
    public Map<String, DuplicateGroup> duplicateGroupMap() { return mDuplicateGroupMap; }

    @VisibleForTesting
    public Set<DuplicateGroup> umiGroups()  { return mDuplicateGroupMap.values().stream().collect(Collectors.toSet()); }

    @VisibleForTesting
    public void processPrimaryFragments(final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList)
    {
        processPrimaryFragments(resolvedFragments, candidateDuplicatesList, Collections.EMPTY_LIST);
    }

    @VisibleForTesting
    public void clearState()
    {
        mFragmentStatus.clear();
        mIncompleteFragments.clear();
        mCandidateDuplicatesMap.clear();
        mDuplicateGroupMap.clear();
    }
}
