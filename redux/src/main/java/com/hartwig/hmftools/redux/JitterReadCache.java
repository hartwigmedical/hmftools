package com.hartwig.hmftools.redux;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.collapseToNonOrientedKeyWithoutCoordinates;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.getFragmentCoordReads;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;
import com.hartwig.hmftools.common.collect.Union;
import com.hartwig.hmftools.common.collect.UnionQuickFind;
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
    private final SortedMultiset<Integer> mInnerCachedReadStarts;
    private final HashMap<String, CollapsedPartition> mCollapsedPartitionLookup = Maps.newHashMap();

    private int mCurrentReadPos;
    private int mLastReadCacheBoundary;

    public JitterReadCache(final ReadCache readCache)
    {
        mReadCache = readCache;
        mMaxSoftClipLength = readCache.maxSoftClipLength();
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

    public void pushSingleRead(final ReadInfo readInfo)
    {
        pushDuplicateGroup(new DuplicateGroup(null, readInfo.read(), readInfo.coordinates()));
    }

    public void pushDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        FragmentCoords coords = duplicateGroup.fragmentCoordinates();
        String collapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);
        mCollapsedPartitionLookup.computeIfAbsent(collapsedKey, key -> new CollapsedPartition());
        mCollapsedPartitionLookup.get(collapsedKey).pushDuplicateGroup(duplicateGroup);
    }

    private FragmentCoordReads popJitterReadCache()
    {
        if(mCollapsedPartitionLookup.isEmpty())
            return null;

        List<DuplicateGroup> groupsToPop = Lists.newArrayList();
        List<String> emptyCollapsedKeys = Lists.newArrayList();
        for(Map.Entry<String, CollapsedPartition> entry : mCollapsedPartitionLookup.entrySet())
        {
            String collapsedKey = entry.getKey();
            CollapsedPartition collapsedPartition = entry.getValue();
            groupsToPop.addAll(collapsedPartition.popReads(mLastReadCacheBoundary));
            if(collapsedPartition.isEmpty())
                emptyCollapsedKeys.add(collapsedKey);
        }

        for(String collapsedKey : emptyCollapsedKeys)
            mCollapsedPartitionLookup.remove(collapsedKey);

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

            return popJitterReadCache();
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
            fragmentCoordReads = popJitterReadCache();

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

        fragmentCoordReads = getFragmentCoordReads(
                mCollapsedPartitionLookup.values().stream().flatMap(CollapsedPartition::duplicateGroupStream));
        mCollapsedPartitionLookup.clear();

        singleReads.addAll(fragmentCoordReads.SingleReads);
        duplicateGroups.addAll(fragmentCoordReads.DuplicateGroups);

        if(singleReads.isEmpty() && duplicateGroups.isEmpty())
            return null;

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    @Override
    public int minCachedReadStart()
    {
        int innerCacheMinReadStart = innerMinCachedReadStart();
        if(innerCacheMinReadStart < 0)
            innerCacheMinReadStart = Integer.MAX_VALUE;

        int outerCacheMinReadStart = mCollapsedPartitionLookup.values()
                .stream()
                .mapToInt(CollapsedPartition::minCachedReadStart)
                .min()
                .orElse(Integer.MAX_VALUE);

        int minReadStart = min(innerCacheMinReadStart, outerCacheMinReadStart);
        return minReadStart == Integer.MAX_VALUE ? -1 : minReadStart;
    }

    @Override
    public int cachedReadCount()
    {
        return mReadCache.cachedReadCount() + mCollapsedPartitionLookup.values()
                .stream()
                .mapToInt(CollapsedPartition::cachedReadCount)
                .sum();
    }

    @Override
    public int cachedFragCoordGroups()
    {
        return mReadCache.cachedFragCoordGroups() + mCollapsedPartitionLookup.values()
                .stream()
                .mapToInt(CollapsedPartition::cachedFragCoordGroups)
                .sum();
    }

    private static List<SAMRecord> allFragmentCoordReads(final FragmentCoordReads fragmentCoordReads)
    {
        if(fragmentCoordReads == null)
            return Lists.newArrayList();

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
        return mCollapsedPartitionLookup.values()
                .stream()
                .flatMap(CollapsedPartition::readNamesStream)
                .collect(Collectors.toCollection(HashMultiset::create));
    }

    private static class CollapsedPartition
    {
        private final static Comparator<Union<Integer, FragmentCoords>> COMPARE_BY_LOWER =
                Comparator.comparingInt(x -> x.hasLeft() ? x.left() : x.right().PositionLower);
        private final static Comparator<Union<Integer, FragmentCoords>> COMPARE_BY_UPPER =
                Comparator.comparingInt(x -> x.hasLeft() ? x.left() : x.right().PositionUpper);

        private final HashMap<FragmentCoords, DuplicateGroup> mDuplicateGroupLookup = Maps.newHashMap();
        private final UnionQuickFind<FragmentCoords> mCoordsMerger = new UnionQuickFind<>();
        private final SortedMap<Integer, TreeSet<Union<Integer, FragmentCoords>>> mLowerToUpper = Maps.newTreeMap();
        private final SortedMap<Integer, TreeSet<Union<Integer, FragmentCoords>>> mUpperToLower = Maps.newTreeMap();
        private final HashMap<FragmentCoords, Integer> mUpperBoundaryLookup = Maps.newHashMap();
        private final TreeMap<Integer, HashMultiset<FragmentCoords>> mUpperBoundaryToFragCoords = Maps.newTreeMap();

        public boolean isEmpty()
        {
            return mDuplicateGroupLookup.isEmpty();
        }

        public int minCachedReadStart()
        {
            return mDuplicateGroupLookup.values().stream()
                    .flatMap(x -> x.reads().stream())
                    .mapToInt(SAMRecord::getAlignmentStart)
                    .min()
                    .orElse(Integer.MAX_VALUE);
        }

        public int cachedReadCount()
        {
            return mDuplicateGroupLookup.values().stream().mapToInt(DuplicateGroup::readCount).sum();
        }

        public int cachedFragCoordGroups()
        {
            return mDuplicateGroupLookup.size();
        }

        public Stream<DuplicateGroup> duplicateGroupStream()
        {
            return mDuplicateGroupLookup.values().stream();
        }

        public Stream<SAMRecord> readStream()
        {
            return duplicateGroupStream().flatMap(x -> x.reads().stream());
        }

        public Stream<String> readNamesStream()
        {
            return readStream().map(SAMRecord::getReadName);
        }

        private void mergeGroups(final FragmentCoords coords1, final FragmentCoords coords2)
        {
            Pair<FragmentCoords, FragmentCoords> fromToPair = mCoordsMerger.merge(coords1, coords2);
            if(fromToPair == null)
                return;

            FragmentCoords fromCoords = fromToPair.getLeft();
            FragmentCoords toCoords = fromToPair.getRight();

            int fromUpperBoundary = mUpperBoundaryLookup.get(fromCoords);
            int toUpperBoundary = mUpperBoundaryLookup.get(toCoords);
            int mergedUpperBoundary = max(fromUpperBoundary, toUpperBoundary);
            mUpperBoundaryLookup.remove(fromCoords);
            mUpperBoundaryLookup.put(toCoords, mergedUpperBoundary);

            mUpperBoundaryToFragCoords.get(fromUpperBoundary).remove(fromCoords);
            mUpperBoundaryToFragCoords.get(toUpperBoundary).remove(toCoords);
            mUpperBoundaryToFragCoords.computeIfAbsent(mergedUpperBoundary, key -> HashMultiset.create());
            mUpperBoundaryToFragCoords.get(mergedUpperBoundary).add(toCoords);

            if(mUpperBoundaryToFragCoords.get(fromUpperBoundary).isEmpty())
                mUpperBoundaryToFragCoords.remove(fromUpperBoundary);

            if(mUpperBoundaryToFragCoords.get(toUpperBoundary).isEmpty())
                mUpperBoundaryToFragCoords.remove(toUpperBoundary);
        }

        public void pushDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            FragmentCoords coords = duplicateGroup.fragmentCoordinates();
            int upperBoundary = coords.readPosition() + SINGLE_END_JITTER_COLLAPSE_DISTANCE;
            mDuplicateGroupLookup.put(coords, duplicateGroup);
            mCoordsMerger.add(coords);
            mUpperBoundaryLookup.put(coords, upperBoundary);
            mUpperBoundaryToFragCoords.computeIfAbsent(upperBoundary, key -> HashMultiset.create());
            mUpperBoundaryToFragCoords.get(upperBoundary).add(coords);

            int posLower = coords.PositionLower;
            int posUpper = coords.PositionUpper;

            mLowerToUpper.computeIfAbsent(posLower, key -> Sets.newTreeSet(COMPARE_BY_UPPER));
            mUpperToLower.computeIfAbsent(posUpper, key -> Sets.newTreeSet(COMPARE_BY_LOWER));

            TreeSet<Union<Integer, FragmentCoords>> uppersFixedLower = mLowerToUpper.get(posLower);
            Union<Integer, FragmentCoords> splitter = Union.createLeft(posUpper);
            Union<Integer, FragmentCoords> floor = uppersFixedLower.floor(splitter);
            Union<Integer, FragmentCoords> ceil = uppersFixedLower.ceiling(splitter);
            if(floor != null && abs(posUpper - floor.right().PositionUpper) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                mergeGroups(coords, floor.right());

            if(ceil != null && abs(posUpper - ceil.right().PositionUpper) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                mergeGroups(coords, ceil.right());

            TreeSet<Union<Integer, FragmentCoords>> lowersFixedUpper = mUpperToLower.get(posUpper);
            splitter = Union.createLeft(posLower);
            floor = lowersFixedUpper.floor(splitter);
            ceil = lowersFixedUpper.ceiling(splitter);
            if(floor != null && abs(posLower - floor.right().PositionLower) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                mergeGroups(coords, floor.right());

            if(ceil != null && abs(posLower - ceil.right().PositionLower) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                mergeGroups(coords, ceil.right());

            Union<Integer, FragmentCoords> coordsForSet = Union.createRight(coords);
            mLowerToUpper.get(posLower).add(coordsForSet);
            mUpperToLower.get(posUpper).add(coordsForSet);
        }

        public Collection<DuplicateGroup> popReads(final int lastReadCacheBoundary)
        {
            List<FragmentCoords> reprCoordsToPop =
                    mUpperBoundaryToFragCoords.headMap(lastReadCacheBoundary).values().stream().flatMap(HashMultiset::stream).toList();
            List<DuplicateGroup> groupsToPop = Lists.newArrayList();
            for(FragmentCoords reprCoords : reprCoordsToPop)
            {
                Collection<FragmentCoords> coordsSet = mCoordsMerger.removeSet(reprCoords);
                for(FragmentCoords coords : coordsSet)
                    groupsToPop.add(mDuplicateGroupLookup.remove(coords));

                mLowerToUpper.get(reprCoords.PositionLower).remove(Union.createLeft(reprCoords.PositionUpper));
                mUpperToLower.get(reprCoords.PositionUpper).remove(Union.createLeft(reprCoords.PositionLower));

                if(mLowerToUpper.get(reprCoords.PositionLower).isEmpty())
                    mLowerToUpper.remove(reprCoords.PositionLower);

                if(mUpperToLower.get(reprCoords.PositionUpper).isEmpty())
                    mUpperToLower.remove(reprCoords.PositionUpper);

                int upperBoundary = mUpperBoundaryLookup.remove(reprCoords);
                mUpperBoundaryToFragCoords.get(upperBoundary).remove(reprCoords);
                if(mUpperBoundaryToFragCoords.get(upperBoundary).isEmpty())
                    mUpperBoundaryToFragCoords.remove(upperBoundary);
            }

            return groupsToPop;
        }
    }
}
