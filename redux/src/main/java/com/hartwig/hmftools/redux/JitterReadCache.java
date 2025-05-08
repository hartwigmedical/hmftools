package com.hartwig.hmftools.redux;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.collapseToNonOrientedKeyWithoutCoordinates;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.getFragmentCoordReads;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;
import com.hartwig.hmftools.common.collect.UnionFind;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMRecord;

public class JitterReadCache implements IReadCache
{
    private final ReadCache mReadCache;
    private final int mMaxSoftClipLength;

    private final Map<String, Set<FragmentCoords>> mAuxFragmentCoordsCache;
    private final Map<FragmentCoords, DuplicateGroup> mAuxDuplicateGroupLookup;

    private final SortedMultiset<Integer> mInnerCachedReadStarts;

    private int mCurrentReadPos;
    private int mLastReadCacheBoundary;

    public JitterReadCache(final ReadCache readCache)
    {
        mReadCache = readCache;
        mMaxSoftClipLength = readCache.maxSoftClipLength();
        mAuxFragmentCoordsCache = Maps.newHashMap();
        mAuxDuplicateGroupLookup = Maps.newHashMap();
        mInnerCachedReadStarts = TreeMultiset.create();

        mCurrentReadPos = 0;
        mLastReadCacheBoundary = readCacheBoundary();
    }

    @Override
    public void processRead(final SAMRecord read)
    {
        mCurrentReadPos = read.getAlignmentStart();
        mInnerCachedReadStarts.add(read.getAlignmentStart());
        mReadCache.processRead(read);
    }

    @Override
    public int currentReadMinPosition()
    {
        return mReadCache.currentReadMinPosition();
    }

    private void pushSingleRead(final ReadInfo readInfo)
    {
        pushDuplicateGroup(new DuplicateGroup(null, readInfo.read(), readInfo.coordinates()));
    }

    private void pushDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        FragmentCoords coords = duplicateGroup.fragmentCoordinates();
        String collapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);
        mAuxFragmentCoordsCache.computeIfAbsent(collapsedKey, x -> Sets.newHashSet());
        mAuxFragmentCoordsCache.get(collapsedKey).add(coords);
        mAuxDuplicateGroupLookup.put(coords, duplicateGroup);
    }

    private int innerMinCachedReadStart()
    {
        if(mInnerCachedReadStarts.isEmpty())
            return -1;

        return mInnerCachedReadStarts.firstEntry().getElement();
    }

    @VisibleForTesting
    public int readCacheBoundary()
    {
        int readCacheMinCachedReadStart = innerMinCachedReadStart();
        if(readCacheMinCachedReadStart < 0)
            return mCurrentReadPos - mMaxSoftClipLength + 1;

        return min(mCurrentReadPos, readCacheMinCachedReadStart) - mMaxSoftClipLength + 1;
    }

    private FragmentCoordReads popAuxCache()
    {
        if(mAuxDuplicateGroupLookup.isEmpty())
            return null;

        List<DuplicateGroup> groupsToPop = Lists.newArrayList();
        List<Pair<String, FragmentCoords>> coordsToDrop = Lists.newArrayList();
        for(String key : mAuxFragmentCoordsCache.keySet())
        {
            List<FragmentCoords> allCoords = Lists.newArrayList();
            UnionFind<FragmentCoords> groupMerger = new UnionFind<>();
            for(FragmentCoords coords : mAuxFragmentCoordsCache.get(key))
            {
                allCoords.add(coords);
                groupMerger.add(coords);
            }

            for(int i = 0; i < allCoords.size() - 1; i++)
            {
                FragmentCoords coords1 = allCoords.get(i);
                for(int j = i + 1; j < allCoords.size(); j++)
                {
                    FragmentCoords coords2 = allCoords.get(j);
                    int lowerDiff = abs(coords1.PositionLower - coords2.PositionLower);
                    int upperDiff = abs(coords1.PositionUpper - coords2.PositionUpper);
                    if(lowerDiff > 0 && upperDiff > 0)
                        continue;

                    if(lowerDiff > SINGLE_END_JITTER_COLLAPSE_DISTANCE || upperDiff > SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                        continue;

                    groupMerger.merge(coords1, coords2);
                }
            }

            for(Set<FragmentCoords> partitions : groupMerger.getPartitions())
            {
                int maxPos = partitions.stream().mapToInt(FragmentCoords::readPosition).max().orElse(-1);
                if(maxPos + SINGLE_END_JITTER_COLLAPSE_DISTANCE >= mLastReadCacheBoundary)
                    continue;

                partitions.forEach(coords ->
                {
                    groupsToPop.add(mAuxDuplicateGroupLookup.get(coords));
                    coordsToDrop.add(Pair.of(key, coords));
                });
            }
        }

        for(Pair<String, FragmentCoords> keyCoordPair : coordsToDrop)
        {
            String key = keyCoordPair.getLeft();
            FragmentCoords coords = keyCoordPair.getRight();
            mAuxFragmentCoordsCache.get(key).remove(coords);
            if(mAuxFragmentCoordsCache.get(key).isEmpty())
                mAuxFragmentCoordsCache.remove(key);

            mAuxDuplicateGroupLookup.remove(coords);
        }

        if(groupsToPop.isEmpty())
            return null;

        return getFragmentCoordReads(groupsToPop);
    }

    @Override
    public FragmentCoordReads popReads()
    {
        boolean tryAuxPop = false;
        FragmentCoordReads fragmentCoordReads = mReadCache.popReads();
        if(fragmentCoordReads != null)
        {
            for(SAMRecord read : allFragmentCoordReads(fragmentCoordReads))
                mInnerCachedReadStarts.remove(read.getAlignmentStart());
        }

        int currentReadCacheBoundary = readCacheBoundary();
        if(mLastReadCacheBoundary != currentReadCacheBoundary)
        {
            mLastReadCacheBoundary = currentReadCacheBoundary;
            tryAuxPop = true;
        }

        if(fragmentCoordReads == null)
        {
            if(!tryAuxPop)
                return null;

            return popAuxCache();
        }

        List<ReadInfo> singleReads = Lists.newArrayList();
        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();
        for(ReadInfo singleRead : fragmentCoordReads.SingleReads)
        {
            if(singleRead.coordinates().PositionUpper == NO_POSITION)
            {
                singleReads.add(singleRead);
                continue;
            }

            if(singleRead.coordinates().SuppReadInfo != null)
            {
                singleReads.add(singleRead);
                continue;
            }

            tryAuxPop = true;
            pushSingleRead(singleRead);
        }

        for(DuplicateGroup duplicateGroup : fragmentCoordReads.DuplicateGroups)
        {
            if(duplicateGroup.fragmentCoordinates().PositionUpper == NO_POSITION)
            {
                duplicateGroups.add(duplicateGroup);
                continue;
            }

            if(duplicateGroup.fragmentCoordinates().SuppReadInfo != null)
            {
                duplicateGroups.add(duplicateGroup);
                continue;
            }

            tryAuxPop = true;
            pushDuplicateGroup(duplicateGroup);
        }

        fragmentCoordReads = null;
        if(tryAuxPop)
            fragmentCoordReads = popAuxCache();

        if(fragmentCoordReads != null)
        {
            singleReads.addAll(fragmentCoordReads.SingleReads);
            duplicateGroups.addAll(fragmentCoordReads.DuplicateGroups);
        }

        if(singleReads.isEmpty() && duplicateGroups.isEmpty())
            return null;

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    @Override
    public FragmentCoordReads evictAll()
    {
        mCurrentReadPos = 0;
        mLastReadCacheBoundary = readCacheBoundary();

        List<ReadInfo> singleReads = Lists.newArrayList();
        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();

        FragmentCoordReads fragmentCoordReads = mReadCache.evictAll();
        mInnerCachedReadStarts.clear();
        if(fragmentCoordReads != null)
        {
            singleReads.addAll(fragmentCoordReads.SingleReads);
            duplicateGroups.addAll(fragmentCoordReads.DuplicateGroups);
        }

        fragmentCoordReads = getFragmentCoordReads(mAuxDuplicateGroupLookup.values());
        mAuxFragmentCoordsCache.clear();
        mAuxDuplicateGroupLookup.clear();
        singleReads.addAll(fragmentCoordReads.SingleReads);
        duplicateGroups.addAll(fragmentCoordReads.DuplicateGroups);

        if(singleReads.isEmpty() && duplicateGroups.isEmpty())
            return null;

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    @Override
    public int minCachedReadStart()
    {
        int minReadStart = innerMinCachedReadStart();
        for(DuplicateGroup group : mAuxDuplicateGroupLookup.values())
        {
            for(SAMRecord read : group.reads())
            {
                int readStart = read.getAlignmentStart();
                if(minReadStart < 0 || readStart < minReadStart)
                    minReadStart = readStart;
            }
        }

        return minReadStart;
    }

    @Override
    public int cachedReadCount()
    {
        return mReadCache.cachedReadCount() + mAuxDuplicateGroupLookup.values().stream().mapToInt(DuplicateGroup::readCount).sum();
    }

    @Override
    public int cachedFragCoordGroups()
    {
        return mReadCache.cachedFragCoordGroups() + mAuxDuplicateGroupLookup.size();
    }

    private static List<SAMRecord> allFragmentCoordReads(final FragmentCoordReads fragmentCoordReads)
    {
        List<SAMRecord> reads = Lists.newArrayList();
        for(DuplicateGroup group : fragmentCoordReads.DuplicateGroups)
            reads.addAll(group.allReads());

        for(ReadInfo readInfo : fragmentCoordReads.SingleReads)
            reads.add(readInfo.read());

        return reads;
    }

    @VisibleForTesting
    public Multiset<String> auxCacheReadNames()
    {
        return mAuxDuplicateGroupLookup.values()
                .stream()
                .flatMap(x -> x.reads().stream())
                .map(SAMRecord::getReadName)
                .collect(Collectors.toCollection(HashMultiset::create));
    }
}
