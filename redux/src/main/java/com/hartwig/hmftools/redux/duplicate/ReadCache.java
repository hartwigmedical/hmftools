package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.floorDiv;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.isSbx;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.redux.common.ReadInfo;

import htsjdk.samtools.SAMRecord;

public class ReadCache implements IReadCache
{
    private final int mGroupSize;
    private final int mMaxSoftClipLength;
    private final boolean mUseFragmentOrientation;
    private final SbxDuplicateCollapser mSbxDuplicateCollapser;

    private int mCurrentReadMinPosition;
    private String mCurrentChromosome;
    private final List<ReadPositionGroup> mPositionGroups;

    private int mLastPopPositionCheck;

    private int mCheckSizeReadCount;
    private int mLastCacheReadCount;

    public static final int DEFAULT_GROUP_SIZE = 200; // larger than the maximum soft-clip length for 151-base reads
    public static final int DEFAULT_MAX_SOFT_CLIP = 150;

    private static final int POP_DISTANCE_CHECK = 100;

    private static final int CHECK_CACHE_READ_COUNT = 10000;
    private static final int LOG_READ_COUNT_THRESHOLD = 100000;

    public ReadCache(
            int groupSize, int maxSoftClipLength, boolean useFragmentOrientation, final DuplicatesConfig duplicatesConfig)
    {
        mGroupSize = groupSize;
        mMaxSoftClipLength = maxSoftClipLength;
        mUseFragmentOrientation = useFragmentOrientation;

        mSbxDuplicateCollapser = isSbx() && duplicatesConfig.SbxMaxDuplicateDistance > 0 ?
                new SbxDuplicateCollapser(duplicatesConfig.SbxMaxDuplicateDistance) : null;

        mPositionGroups = Lists.newArrayList();
        mCurrentReadMinPosition = 0;
        mCurrentChromosome = "";
        mLastPopPositionCheck = 0;
        mLastCacheReadCount = 0;
        mCheckSizeReadCount = 0;
    }

    @VisibleForTesting
    public ReadCache(int groupSize, int maxSoftClipLength, boolean useFragmentOrientation)
    {
        this(groupSize, maxSoftClipLength, useFragmentOrientation, new DuplicatesConfig(0));
    }

    @Override
    public void processRead(final SAMRecord read)
    {
        int readLeftSoftClip = leftSoftClipLength(read);

        if(readLeftSoftClip > mMaxSoftClipLength)
        {
            RD_LOGGER.warn("read({}) leftClip({}) exceeds maxSoftClipLength({})",
                    readToString(read), readLeftSoftClip, mMaxSoftClipLength);
        }

        FragmentCoords fragmentCoords = FragmentCoords.fromRead(read, mUseFragmentOrientation);

        String chromosome = read.getReferenceName();
        ReadPositionGroup group = getOrCreateGroup(chromosome, fragmentCoords);

        group.addRead(fragmentCoords, read);

        mCurrentReadMinPosition = read.getAlignmentStart();

        if(!mCurrentChromosome.equals(chromosome))
        {
            mCurrentChromosome = chromosome;
            mLastPopPositionCheck = 0;
        }

        checkCacheSize();
    }

    public int maxSoftClipLength() { return mMaxSoftClipLength; }

    @Override
    public int currentReadMinPosition() { return mCurrentReadMinPosition; }

    @Override
    public FragmentCoordReads popReads()
    {
        return popPassedReads(true);
    }

    @Override
    public FragmentCoordReads evictAll()
    {
        return popPassedReads(false);
    }

    private FragmentCoordReads popPassedReads(boolean checkCurrentPosition)
    {
        // process and remove any fragment coords which:
        // 1. Uses the lower read position (ie unclipped 5' position from a read start)
        //  - if the fragment coords lower position <= max soft-clip length from current read start position
        // 2. Uses the upper read position (ie unclipped read end position)
        //  - if the fragment coords upper position < current read start position
        //
        // then remove any read groups which are then empty

        if(checkCurrentPosition)
        {
            // process if has an earlier chromosome or if the current position has moved past the last check position
            if(!mPositionGroups.isEmpty() && mPositionGroups.get(0).Chromosome.equals(mCurrentChromosome))
            {
                if(mCurrentReadMinPosition < mLastPopPositionCheck + POP_DISTANCE_CHECK)
                    return null;
            }
        }

        mLastPopPositionCheck = (int)(floor(mCurrentReadMinPosition/POP_DISTANCE_CHECK) * POP_DISTANCE_CHECK);

        List<DuplicateGroup> duplicateGroups = null;
        List<ReadInfo> singleReads = null;

        int popFragCoordLowerPosition = mCurrentReadMinPosition - mMaxSoftClipLength;
        int popFragCoordUpperPosition = mCurrentReadMinPosition - 1;

        int groupIndex = 0;
        while(groupIndex < mPositionGroups.size())
        {
            ReadPositionGroup group = mPositionGroups.get(groupIndex);

            if(checkCurrentPosition)
            {
                if(group.Chromosome.equals(mCurrentChromosome) && mCurrentReadMinPosition < group.PositionStart)
                    break;
            }

            boolean takeAllFragments = !group.Chromosome.equals(mCurrentChromosome) || group.PositionEnd < popFragCoordLowerPosition;

            Set<FragmentCoords> processedCoords = Sets.newHashSet();

            for(Map.Entry<FragmentCoords,List<SAMRecord>> entry : group.FragCoordsMap.entrySet())
            {
                FragmentCoords fragCoords = entry.getKey();

                if(checkCurrentPosition && !takeAllFragments && group.Chromosome.equals(mCurrentChromosome))
                {
                    int readPosition = fragCoords.readPosition();
                    Orientation readOrientation = fragCoords.readOrientation();

                    if(readOrientation.isForward())
                    {
                        if(readPosition > popFragCoordLowerPosition)
                            continue;
                    }
                    else
                    {
                        if(readPosition > popFragCoordUpperPosition)
                            continue;
                    }
                }

                processedCoords.add(fragCoords);

                List<SAMRecord> reads = entry.getValue();

                if(reads.size() > 1)
                {
                    if(duplicateGroups == null)
                        duplicateGroups = Lists.newArrayList();

                    duplicateGroups.add(new DuplicateGroup(reads, fragCoords));
                }
                else
                {
                    if(singleReads == null)
                        singleReads = Lists.newArrayList();

                    singleReads.add(new ReadInfo(reads.get(0), fragCoords));
                }
            }

            if(!processedCoords.isEmpty())
            {
                processedCoords.forEach(x -> group.FragCoordsMap.remove(x));

                if(group.FragCoordsMap.isEmpty())
                {
                    mPositionGroups.remove(groupIndex);
                    continue;
                }
            }

            ++groupIndex;
        }

        if(duplicateGroups == null && singleReads == null)
            return null;

        if(mSbxDuplicateCollapser == null)
            return new FragmentCoordReads(duplicateGroups, singleReads);
        else
            return mSbxDuplicateCollapser.collapseGroups(duplicateGroups, singleReads);
    }

    @Override
    public int minCachedReadStart()
    {
        return mPositionGroups.stream().mapToInt(x -> x.minReadPosition()).min().orElse(-1);
    }

    private ReadPositionGroup getOrCreateGroup(final String chromosome, final FragmentCoords fragmentCoords)
    {
        // find an existing group by the applicable fragment coordination position, otherwise make a new one
        int groupIndex = 0;

        int fragmentPosition = fragmentCoords.readPosition();

        if(!mPositionGroups.isEmpty())
        {
            for(; groupIndex < mPositionGroups.size(); ++groupIndex)
            {
                ReadPositionGroup group = mPositionGroups.get(groupIndex);

                if(!chromosome.equals(group.Chromosome))
                    continue;

                if(positionWithin(fragmentPosition, group.PositionStart, group.PositionEnd))
                    return group;

                if(fragmentPosition < group.PositionStart)
                    break;
            }
        }

        int remainder = fragmentPosition % mGroupSize;

        int groupPosStart, groupPosEnd;

        if(remainder == 0)
        {
            groupPosEnd = mGroupSize * floorDiv(fragmentPosition, mGroupSize);
            groupPosStart = groupPosEnd - mGroupSize + 1;
        }
        else
        {
            groupPosStart = mGroupSize * floorDiv(fragmentPosition, mGroupSize) + 1;
            groupPosEnd = groupPosStart + mGroupSize - 1;
        }

        ReadPositionGroup group = new ReadPositionGroup(chromosome, groupPosStart, groupPosEnd);
        mPositionGroups.add(groupIndex, group);
        return group;
    }

    private class ReadPositionGroup
    {
        public final String Chromosome;
        public final int PositionStart;
        public final int PositionEnd;

        public final Map<FragmentCoords,List<SAMRecord>> FragCoordsMap;

        public ReadPositionGroup(final String chromosome, final int positionStart, final int positionEnd)
        {
            Chromosome = chromosome;
            PositionStart = positionStart;
            PositionEnd = positionEnd;
            FragCoordsMap = Maps.newHashMap();
        }

        public void addRead(final FragmentCoords fragCoords, final SAMRecord read)
        {
            List<SAMRecord> reads = FragCoordsMap.get(fragCoords);

            if(reads == null)
            {
                reads = Lists.newArrayList();
                FragCoordsMap.put(fragCoords, reads);
            }

            reads.add(read);
        }

        public int readCount() { return FragCoordsMap.values().stream().mapToInt(x -> x.size()).sum(); }

        public int minReadPosition()
        {
            int minReadPosition = -1;

            for(List<SAMRecord> reads : FragCoordsMap.values())
            {
                int minReadPos = reads.stream().mapToInt(x -> x.getAlignmentStart()).min().orElse(-1);
                if(minReadPosition < 0 || minReadPos < minReadPosition)
                    minReadPosition = minReadPos;
            }

            return minReadPosition;
        }

        public String toString()
        {
            return format("range(%s:%d-%d) frags(%d) reads(%d)",
                    Chromosome, PositionStart, PositionEnd, FragCoordsMap.size(), readCount());
        }
    }

    public String toString()
    {
        return format("minPos(%s:%d) groups(%d) frags(%d) reads(%d)",
                mCurrentChromosome, mCurrentReadMinPosition, mPositionGroups.size(), cachedFragCoordGroups(), cachedReadCount());
    }

    private void checkCacheSize()
    {
        // only check and log cache size every X reads
        ++mCheckSizeReadCount;

        if(mCheckSizeReadCount < CHECK_CACHE_READ_COUNT)
            return;

        mCheckSizeReadCount = 0;

        int newReadCount = cachedReadCount();

        if(newReadCount < LOG_READ_COUNT_THRESHOLD)
            return;

        if(abs(newReadCount - mLastCacheReadCount) < LOG_READ_COUNT_THRESHOLD)
            return;

        RD_LOGGER.debug("read cache({}) above threshold", toString());
        mLastCacheReadCount = newReadCount;
    }

    @Override
    @VisibleForTesting
    public int cachedReadCount() { return mPositionGroups.stream().mapToInt(x -> x.readCount()).sum(); }

    @Override
    public int cachedFragCoordGroups() { return mPositionGroups.stream().mapToInt(x -> x.FragCoordsMap.size()).sum(); }
    public int cachedReadGroups() { return mPositionGroups.size(); }

    public void clear()
    {
        mPositionGroups.clear();
        mCurrentChromosome = "";
        mCurrentReadMinPosition = 0;
        mLastPopPositionCheck = 0;
        mLastCacheReadCount = 0;
        mCheckSizeReadCount = 0;
    }
}
