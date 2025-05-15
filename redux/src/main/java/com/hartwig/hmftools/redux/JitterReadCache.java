package com.hartwig.hmftools.redux;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;
import com.hartwig.hmftools.redux.common.AuxReadCache_0;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.IAuxReadCache;
import com.hartwig.hmftools.redux.common.ReadInfo;

import htsjdk.samtools.SAMRecord;

public class JitterReadCache implements IReadCache
{
    private final ReadCache mReadCache;
    private final int mMaxSoftClipLength;
    private final IAuxReadCache mAuxReadCache = new AuxReadCache_0();
    private final SortedMultiset<Integer> mInnerCachedReadStarts;

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

            return mAuxReadCache.popReads(mLastReadCacheBoundary);
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
            mAuxReadCache.pushSingleRead(singleRead);
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
            mAuxReadCache.pushDuplicateGroup(duplicateGroup);
        }

        fragmentCoordReads = null;
        if(tryAuxPop)
            fragmentCoordReads = mAuxReadCache.popReads(mLastReadCacheBoundary);

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

        fragmentCoordReads = mAuxReadCache.evictAll();

        singleReads.addAll(fragmentCoordReads.SingleReads);
        duplicateGroups.addAll(fragmentCoordReads.DuplicateGroups);

        if(singleReads.isEmpty() && duplicateGroups.isEmpty())
            return null;

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    @Override
    public int minCachedReadStart()
    {
        int innerMinReadStart = innerMinCachedReadStart();
        int auxMinReadStart = mAuxReadCache.minCachedReadStart();
        if(innerMinReadStart < 0 && auxMinReadStart < 0)
            return -1;

        if(innerMinReadStart < 0)
            innerMinReadStart = Integer.MAX_VALUE;

        if(auxMinReadStart < 0)
            auxMinReadStart = Integer.MAX_VALUE;

        return min(innerMinReadStart, auxMinReadStart);
    }

    @Override
    public int cachedReadCount()
    {
        return mReadCache.cachedReadCount() + mAuxReadCache.cachedReadCount();
    }

    @Override
    public int cachedFragCoordGroups() { return mReadCache.cachedFragCoordGroups() + mAuxReadCache.cachedFragCoordGroups(); }

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
    public Multiset<String> auxCacheReadNames() { return mAuxReadCache.readNames(); }
}
