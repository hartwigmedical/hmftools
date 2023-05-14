package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_CONSENSUS_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.NANOS_IN_SECOND;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.overlapsExcludedRegion;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;
import static com.hartwig.hmftools.markdups.common.ReadMatch.NO_READ_MATCH;
import static com.hartwig.hmftools.markdups.common.ResolvedFragmentState.fragmentState;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.markdups.RecordWriter;
import com.hartwig.hmftools.markdups.umi.ConsensusReads;
import com.hartwig.hmftools.markdups.umi.UmiConfig;
import com.hartwig.hmftools.markdups.umi.UmiGroup;

import htsjdk.samtools.SAMRecord;

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
    private long mLockAcquireTime;
    private boolean mPerfChecks;

    private Set<UmiGroup> mUpdatedUmiGroups;
    private Set<CandidateDuplicates> mUpdatedCandidateDuplicates;
    private ChrBaseRegion mExcludedRegion;

    private static final int LOG_CACHE_COUNT = 10000;

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
        mExcludedRegion = null;

        mLock = new ReentrantLock();
        mLockAcquireTime = 0;
        mPerfChecks = false;
    }

    public String partitionStr() { return mChrPartition; }
    public Statistics statistics() { return mDuplicateGroups.statistics(); }

    public void togglePerfChecks() { mPerfChecks = true; }
    public double totalLockTime() { return mLockAcquireTime / NANOS_IN_SECOND; }

    public synchronized void setExcludedRegion(final ChrBaseRegion excludedRegion) { mExcludedRegion = excludedRegion; }

    public void processPrimaryFragments(
            final List<Fragment> resolvedFragments, final List<CandidateDuplicates> candidateDuplicatesList, final List<UmiGroup> umiGroups)
    {
        // gather any cached mate reads, attempt to resolve any candidate duplicates and feed back the resultant set of resolved fragments
        try
        {
            acquireLock();

            // UMIs are filtered here since the fragments themselves don't need to collect incomplete reads nor set resolved status
            resolvedFragments.stream().filter(x -> x.umi() == null).forEach(x -> processResolvedFragment(x));

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

        UmiGroup umiGroup = mUmiGroups.get(read.getReadName());

        if(umiGroup != null)
        {
            umiGroup.addRead(read);
            mUpdatedUmiGroups.add(umiGroup);
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
            MD_LOGGER.error("umiGroup({}) has no read IDs: {}", umiGroup.id(), umiGroup.toString());
            return;
        }

        groupReadIds.forEach(x -> mUmiGroups.remove(x));
    }

    private void checkResolveCandidateDuplicates(final CandidateDuplicates candidateDuplicates)
    {
        if(!candidateDuplicates.allFragmentsReady())
            return;

        List<List<Fragment>> duplicateGroups = candidateDuplicates.finaliseFragmentStatus();

        boolean inExcludedRegion = mExcludedRegion != null && duplicateGroups.stream()
                .anyMatch(x -> x.stream().anyMatch(y -> y.reads().stream().anyMatch(z -> overlapsExcludedRegion(mExcludedRegion, z))));

        if(umiEnabled() && !inExcludedRegion)
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
            mDuplicateGroups.processDuplicateGroups(duplicateGroups, inExcludedRegion);
        }

        for(Fragment fragment : candidateDuplicates.fragments())
        {
            mIncompleteFragments.remove(fragment.id());

            // store this new resolved state if more reads are expected for the fragment
            if(fragment.allReadsPresent())
                continue;

            if(fragment.umi() != null) // cached with the UMI group
                continue;

            ResolvedFragmentState resolvedState = fragmentState(fragment);
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

    public int writeRemainingReads(final RecordWriter recordWriter, final ConsensusReads consensusReads, boolean logCachedReads)
    {
        if(mUmiGroups.isEmpty() && mIncompleteFragments.isEmpty() && mFragmentStatus.isEmpty())
            return 0;

        if(logCachedReads && MD_LOGGER.isDebugEnabled())
        {
            // not under lock since called only when all partitions are complete
            MD_LOGGER.debug("partition({}) final state: {}", mChrPartition, cacheCountsStr());
        }

        // a clean-up routine for cached fragments
        Set<UmiGroup> processedUmiGroups = Sets.newHashSet();

        int cachedReadCount = 0;

        for(UmiGroup umiGroup : mUmiGroups.values())
        {
            if(processedUmiGroups.contains(umiGroup))
                continue;

            processedUmiGroups.add(umiGroup);

            for(String readId : umiGroup.getReadIds())
            {
                Fragment incompleteFragment = mIncompleteFragments.get(readId);

                if(incompleteFragment != null)
                {
                    mIncompleteFragments.remove(readId);
                    incompleteFragment.reads().forEach(x -> umiGroup.addRead(x));
                }
            }

            int cachedUmiReads = umiGroup.cachedReadCount();

            if(cachedUmiReads == 0)
                continue;

            cachedReadCount += cachedUmiReads;

            List<SAMRecord> completeReads = umiGroup.popCompletedReads(consensusReads, true);
            recordWriter.writeUmiReads(umiGroup, completeReads);

            if(logCachedReads)
            {
                MD_LOGGER.debug("writing {} cached reads for umi group({}) coords({})",
                        cachedUmiReads, umiGroup.toString(), umiGroup.coordinatesKey());

                for(SAMRecord read : completeReads)
                {
                    if(read.getSupplementaryAlignmentFlag() || read.hasAttribute(UMI_CONSENSUS_ATTRIBUTE))
                        continue;

                    MD_LOGGER.debug("writing umi read: {}", readToString(read));
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
                    MD_LOGGER.debug("writing incomplete read: {} status({})", readToString(read), fragment.status());
                }
            }
        }

        if(logCachedReads && !mFragmentStatus.isEmpty())
        {
            for(Map.Entry<String,ResolvedFragmentState> entry : mFragmentStatus.entrySet())
            {
                MD_LOGGER.debug("cached resolved status: {} : {}", entry.getKey(), entry.getValue());
            }
        }

        mFragmentStatus.clear();
        mIncompleteFragments.clear();
        mCandidateDuplicatesMap.clear();
        mUmiGroups.clear();

        return cachedReadCount;
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

        long umiReads = 0;
        Set<UmiGroup> uniqueGroups = Sets.newHashSet();

        if(mPerfChecks)
        {
            mUmiGroups.values().forEach(x -> uniqueGroups.add(x));
            umiReads = uniqueGroups.stream().mapToInt(x -> x.cachedReadCount()).sum();
        }

        return format("incomplete(%d supp=%d) resolved(%d supp=%d mate=%d) umi(groups=%d frags=%s reads=%d) candidateGroups(%d max=%d)",
                mIncompleteFragments.size(), incompleteSupp, mFragmentStatus.size(), resolvedNoSupp, resolvedNoMate,
                uniqueGroups.size(), mUmiGroups.size(), umiReads, mCandidateDuplicatesMap.size(), maxCandidateGroup);
    }

    public void logCacheCounts()
    {
        try
        {
            acquireLock();

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
            acquireLock();

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
                mChrPartition, mFragmentStatus.size(), mIncompleteFragments.size(), mCandidateDuplicatesMap.size(), mUmiGroups.size());
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
