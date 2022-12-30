package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.chromosomeIndicator;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

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
            final String chrPartition, final List<Fragment> resolvedFragments, final List<Fragment> supplementaries)
    {
        // cache entries which are expecting further reads, otherwise combine reads from existing fragments and set resolved status

        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        for(Fragment fragment : resolvedFragments)
        {
            handleResolvedFragment(fragment, chrPartition);
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

        if(mMultiChromosomes && !fragment.hasRemotePartitions())
            return;

        String fragmentChromosome = chromosomeIndicator(fragment.reads().get(0).getReferenceName());

        List<String> otherPartitions = mMultiChromosomes ? fragment.remotePartitions() : Lists.newArrayList(chrPartition);

        for(String remotePartition : otherPartitions)
        {
            // local GC stores and applies same chromosome partitions only and vice versa
            if(mSingleChromosome && !remotePartition.startsWith(fragmentChromosome)) // only applicable for remote GC
                continue;
            else if(mMultiChromosomes && remotePartition.startsWith(fragmentChromosome)) // already applied in local GC
                continue;

            boolean isProcessed = mProcessedPartitions.contains(remotePartition);

            boolean isCurrentPartition = mSingleChromosome && remotePartition.equals(chrPartition);

            if(isProcessed || isCurrentPartition)
            {
                setResolvedStatus(remotePartition, fragment);
            }

            if(!isProcessed || isCurrentPartition)
            {
                // store status against the remote partition, not the current one
                storeResolvedStatus(remotePartition, fragment);
            }
        }
    }

    private void setResolvedStatus(final String chrPartition, final Fragment resolvedFragment)
    {
        // when a fragment is resolved and there are unclear fragments for the same read ID, then write them with the resolved status
        PartitionCache partitionCache = mPartitionCacheMap.get(chrPartition);

        if(partitionCache == null)
            return;

        Fragment fragment = partitionCache.Supplementaries.get(resolvedFragment.id());

        if(fragment == null)
            return;

        fragment.setStatus(resolvedFragment.status());
        mRecordWriter.writeFragment(fragment);
        partitionCache.Supplementaries.remove(resolvedFragment.id());
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
        String remotePartition = mMultiChromosomes ? supplementary.remotePartitions().get(0) : partitionCache.ChrPartition;

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
            Fragment fragment = partitionCache.Supplementaries.get(supplementary.id());

            if(fragment != null)
                fragment.merge(supplementary);
            else
                partitionCache.Supplementaries.put(supplementary.id(), supplementary);
        }
    }

    public void localPartitionComplete(final String chrPartition)
    {
        partitionComplete(chrPartition);
    }

    private void partitionComplete(final String chrPartition)
    {
        mProcessedPartitions.add(chrPartition);

        // write unmatched local supplementaries
        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        // remove fragment status since no longer required
        partitionCache.FragmentStatus.clear();

        if(mSingleChromosome)
        {
            for(Fragment fragment : partitionCache.Supplementaries.values())
            {
                if(fragment.status() == SUPPLEMENTARY)
                {
                    mRecordWriter.writeFragment(fragment);
                    BM_LOGGER.debug("supplementary({}) in local GC unmatched", readToString(fragment.reads().get(0)));
                }
            }

            partitionCache.Supplementaries.clear();
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

    /*
    public void gatherPartitionFragments(
            final String chrPartition, final BaseRegion partitionRegion, final List<Fragment> resolvedFragments)
    {
        if(mMultiChromosomes)
        {
            BM_LOGGER.error("remote GC called local function for partition({})", chrPartition);
            return;
        }

        PartitionCache partitionCache = getOrCreatePartitionCache(chrPartition);

        for(Fragment fragment : partitionCache.Supplementaries.values())
        {
            if(fragment.status() == SUPPLEMENTARY)
            {
                mRecordWriter.writeFragment(fragment);
                BM_LOGGER.debug("supplementary({}) in local GC unmatched", readToString(fragment.reads().get(0)));
            }
        }

        partitionCache.Supplementaries.clear();

        // partitionCache.clear(); // keep resolved status for subsequent reads

        mProcessedPartitions.add(chrPartition);
    }
    */

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
                    mRecordWriter.writeResolvedReadData(entry.getKey(), entry.getValue(), partitionCache.ChrPartition);
                }
            }
            else
            {
                for(Fragment fragment : partitionCache.Supplementaries.values())
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