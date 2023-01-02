package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.chromosomeIndicator;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class GroupCombiner
{
    private final boolean mSingleChromosome;
    private final boolean mMultiChromosomes;
    private final boolean mWriteCandidates;
    private final RecordWriter mRecordWriter;

    private final Map<String, PartitionCache> mPartitionCacheMap;
    private final Set<String> mProcessedPartitions;

    public GroupCombiner(final RecordWriter recordWriter, boolean singleChromosome, boolean writeCandidates)
    {
        mRecordWriter = recordWriter;
        mPartitionCacheMap = Maps.newHashMap();
        mProcessedPartitions = Sets.newHashSet();
        mSingleChromosome = singleChromosome;
        mMultiChromosomes = !mSingleChromosome;
        mWriteCandidates = writeCandidates;
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
            final String chrPartition, final List<Fragment> resolvedFragments,
            final List<CandidateDuplicates> candidateDuplicatesList,
            final List<Fragment> supplementaries)
    {
        // cache entries which are expecting further reads, otherwise combine reads from existing fragments and set resolved status

        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        for(Fragment fragment : resolvedFragments)
        {
            handleResolvedFragment(fragment, chrPartition);
        }

        for(CandidateDuplicates candidateDuplicates : candidateDuplicatesList)
        {
            // handleCandidateDuplicates(candidateDuplicates, partitionCache);
        }

        for(Fragment fragment : supplementaries)
        {
            handleSupplementary(fragment, partitionCache);
        }

        if(mMultiChromosomes)
            partitionComplete(chrPartition);
    }

    public FragmentStatus findFragmentStatus(final String chrPartition, final String readId)
    {
        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);
        FragmentStatus status = partitionCache.FragmentStatus.get(readId);
        return status != null ? status : FragmentStatus.UNSET;
    }

    private void handleResolvedFragment(final Fragment fragment, final String chrPartition)
    {
        if(fragment.allReadsPresent())
            return;

        if(!fragment.hasRemotePartitions())
            return;

        String fragmentChromosome = chromosomeIndicator(fragment.reads().get(0).getReferenceName());

        for(String remotePartition : fragment.remotePartitions())
        {
            // ignore the local partition which may also have a supplementary for this fragment
            if(mMultiChromosomes && remotePartition.startsWith(fragmentChromosome))
                continue;

            if(mProcessedPartitions.contains(remotePartition))
            {
                setResolvedStatus(remotePartition, fragment);
            }
            else
            {
                // store status against the remote partition, not the current one
                storeResolvedStatus(remotePartition, fragment);
            }
        }
    }

    public void localResolvedFragment(final String chrPartition, final Fragment fragment)
    {
        setResolvedStatus(chrPartition, fragment);
        storeResolvedStatus(chrPartition, fragment);
    }

    private void setResolvedStatus(final String chrPartition, final Fragment resolvedFragment)
    {
        // when a fragment is resolved and there are unclear fragments for the same read ID, then write them with the resolved status
        PartitionCache partitionCache = mPartitionCacheMap.get(chrPartition);

        if(partitionCache == null)
            return;

        Fragment fragment = partitionCache.IncompleteFragments.get(resolvedFragment.id());

        if(fragment == null)
            return;

        fragment.setStatus(resolvedFragment.status());
        mRecordWriter.writeFragment(fragment);
        partitionCache.IncompleteFragments.remove(resolvedFragment.id());
    }

    private void storeResolvedStatus(final String chrPartition, final Fragment resolvedFragment)
    {
        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);
        partitionCache.FragmentStatus.put(resolvedFragment.id(), resolvedFragment.status());
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

    public void localSupplementary(final String chrPartition, final Fragment supplementary)
    {
        // not synchronised since only called locally
        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);
        handleSupplementary(supplementary, partitionCache);
    }

    private void handleSupplementary(final Fragment supplementary, final PartitionCache partitionCache)
    {
        String remotePartition = mMultiChromosomes ? supplementary.remotePartitions().get(0) : partitionCache.mChrPartition;

        boolean processedPartition = mProcessedPartitions.contains(remotePartition);

        boolean storeFragment = !processedPartition || mSingleChromosome;

        if(processedPartition || mSingleChromosome)
        {
            // check for a resolved status (always stored against the remote partition)
            if(checkResolvedStatus(supplementary, partitionCache))
                return;

            // no read group came through on this remote partition - shouldn't happen, but store the read to be written at the end
            storeFragment = true;
        }

        if(storeFragment)
        {
            if(mWriteCandidates)
            {
                mRecordWriter.writeCachedFragment(supplementary);
                return;
            }

            // store supplementary for when this partition is processed
            Fragment fragment = partitionCache.IncompleteFragments.get(supplementary.id());

            if(fragment != null)
                supplementary.reads().forEach(x -> fragment.addRead(x));
            else
                partitionCache.IncompleteFragments.put(supplementary.id(), supplementary);
        }
    }

    public void localPartitionComplete(final String chrPartition)
    {
        partitionComplete(chrPartition);
    }

    private void partitionComplete(final String chrPartition)
    {
        mProcessedPartitions.add(chrPartition);

        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        // remove fragment status since no longer required
        partitionCache.FragmentStatus.clear();

        if(mSingleChromosome)
        {
            // write unmatched local supplementaries
            for(Fragment fragment : partitionCache.IncompleteFragments.values())
            {
                if(fragment.status() == SUPPLEMENTARY)
                {
                    mRecordWriter.writeFragment(fragment);
                    BM_LOGGER.debug("supplementary({}) in local GC unmatched", readToString(fragment.reads().get(0)));
                }
            }

            partitionCache.IncompleteFragments.clear();
        }
        else
        {
            if(mMultiChromosomes && (mProcessedPartitions.size() % 100) == 0)
            {
                int cachedSupplementaries = mPartitionCacheMap.values().stream().mapToInt(x -> x.incompleteFragments()).sum();
                int cachedResolved = mPartitionCacheMap.values().stream().mapToInt(x -> x.resolvedFragments()).sum();

                BM_LOGGER.info("remote group cache partitions processed({}) cached supp({}) resolved({})",
                        mProcessedPartitions.size(), cachedSupplementaries, cachedResolved);
            }
        }
    }

    public void handleRemaining()
    {
        int cachedSupplementaries = mPartitionCacheMap.values().stream().mapToInt(x -> x.incompleteFragments()).sum();
        int cachedResolved = mPartitionCacheMap.values().stream().mapToInt(x -> x.resolvedFragments()).sum();

        BM_LOGGER.info("final group cache partitions({}) supps({}) resolved({})",
                mProcessedPartitions.size(), cachedSupplementaries, cachedResolved);

        for(PartitionCache partitionCache : mPartitionCacheMap.values())
        {
            if(mWriteCandidates)
            {
                for(Map.Entry<String,FragmentStatus> entry : partitionCache.FragmentStatus.entrySet())
                {
                    mRecordWriter.writeResolvedReadData(entry.getKey(), entry.getValue(), partitionCache.mChrPartition);
                }
            }
            else
            {
                for(Fragment fragment : partitionCache.IncompleteFragments.values())
                {
                    mRecordWriter.writeFragment(fragment);
                }
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