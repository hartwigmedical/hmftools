package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.collapseToNonOrientedKeyWithoutCoordinates;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.getFragmentCoordReads;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.collect.UnionFind;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMRecord;

// TODO: remove
public class AuxReadCache_0 implements IAuxReadCache
{
    private final Map<String, Set<FragmentCoords>> mAuxFragmentCoordsCache = Maps.newHashMap();
    private final Map<FragmentCoords, DuplicateGroup> mAuxDuplicateGroupLookup = Maps.newHashMap();

    @Override
    public int minCachedReadStart()
    {
        int minReadStart = -1;
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
        return mAuxDuplicateGroupLookup.values().stream().mapToInt(DuplicateGroup::readCount).sum();
    }

    @Override
    public int cachedFragCoordGroups() { return mAuxDuplicateGroupLookup.size(); }

    @Override
    public Multiset<String> readNames()
    {
        return mAuxDuplicateGroupLookup.values().stream().flatMap(x -> x.reads().stream()).map(SAMRecord::getReadName).collect(Collectors.toCollection(HashMultiset::create));
    }

    @Override
    public void pushSingleRead(final ReadInfo readInfo)
    {
        pushDuplicateGroup(new DuplicateGroup(null, readInfo.read(), readInfo.coordinates()));
    }

    @Override
    public void pushDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        FragmentCoords coords = duplicateGroup.fragmentCoordinates();
        String collapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);
        mAuxFragmentCoordsCache.computeIfAbsent(collapsedKey, x -> Sets.newHashSet());
        mAuxFragmentCoordsCache.get(collapsedKey).add(coords);
        mAuxDuplicateGroupLookup.put(coords, duplicateGroup);
    }

    @Override
    public FragmentCoordReads popReads(final int lastReadCacheBoundary)
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
                if(maxPos + SINGLE_END_JITTER_COLLAPSE_DISTANCE >= lastReadCacheBoundary)
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
    public FragmentCoordReads evictAll()
    {
        FragmentCoordReads fragmentCoordReads = getFragmentCoordReads(mAuxDuplicateGroupLookup.values());
        mAuxFragmentCoordsCache.clear();
        mAuxDuplicateGroupLookup.clear();
        return fragmentCoordReads;
    }
}
