package com.hartwig.hmftools.redux;

import static java.lang.Math.floor;
import static java.lang.Math.sin;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.common.ReadInfo;

import htsjdk.samtools.SAMRecord;

public class ReadCache
{
    private final int mGroupSize;
    private int mCurrentReadMinPosition;
    private String mCurrentChromosome;
    private final List<ReadPositionGroup> mPositionGroups;

    public static final int DEFAULT_GROUP_SIZE = 200; // larger than the maximum soft-clip length for 151-base reads

    public ReadCache(int groupSize)
    {
        mGroupSize = groupSize;
        mPositionGroups = Lists.newArrayList();
        mCurrentReadMinPosition = 0;
        mCurrentChromosome = "";
    }

    public void processRead(final SAMRecord read)
    {
        FragmentCoords fragmentCoords = FragmentCoords.fromRead(read);

        ReadPositionGroup group = getOrCreateGroup(fragmentCoords);

        group.addRead(fragmentCoords, read);

        mCurrentReadMinPosition = read.getAlignmentStart();
        mCurrentChromosome = read.getReferenceName();
    }

    public FragmentCoordReads popReadsByFragmentCoordinates()
    {
        if(mPositionGroups.size() < 3)
            return null;

        if(mCurrentReadMinPosition <= mPositionGroups.get(1).PositionEnd)
            return null;

        List<DuplicateGroup> duplicateGroups = null;
        List<ReadInfo> singleReads = null;

        while(!mPositionGroups.isEmpty())
        {
            ReadPositionGroup group = mPositionGroups.get(0);

            // TODO: test chromosome change: mCurrentChromosome

            if(mCurrentReadMinPosition <= group.PositionEnd + mGroupSize)
                break;

            for(Map.Entry<FragmentCoords,List<SAMRecord>> entry : group.FragCoordsMap.entrySet())
            {
                List<SAMRecord> reads = entry.getValue();

                if(reads.size() > 1)
                {
                    if(duplicateGroups == null)
                        duplicateGroups = Lists.newArrayList();

                    duplicateGroups.add(new DuplicateGroup(reads, entry.getKey()));
                }
                else
                {
                    if(singleReads == null)
                        singleReads = Lists.newArrayList();

                    singleReads.add(new ReadInfo(reads.get(0), entry.getKey()));
                }
            }

            mPositionGroups.remove(0);
        }

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    public List<List<SAMRecord>> processRemaining()
    {
        if(mPositionGroups.isEmpty())
            return null;

        List<List<SAMRecord>> fragCoordReadsList = Lists.newArrayList();

        while(!mPositionGroups.isEmpty())
        {
            ReadPositionGroup group = mPositionGroups.get(0);

            for(List<SAMRecord> reads : group.FragCoordsMap.values())
            {
                fragCoordReadsList.add(reads);
            }

            mPositionGroups.remove(0);
        }

        return fragCoordReadsList;
    }

    private ReadPositionGroup getOrCreateGroup(final FragmentCoords fragmentCoords)
    {
        // select the lowest existing group by position
        int groupIndex = 0;

        int fragmentPosition = fragmentCoords.ReadIsLower ? fragmentCoords.PositionLower : fragmentCoords.PositionUpper;

        if(!mPositionGroups.isEmpty())
        {
            for(; groupIndex < mPositionGroups.size(); ++groupIndex)
            {
                ReadPositionGroup group = mPositionGroups.get(groupIndex);

                if(positionWithin(fragmentPosition, group.PositionStart, group.PositionEnd))
                    return group;

                if(fragmentPosition < group.PositionStart)
                    break;
            }
        }

        int groupPosStart = mGroupSize * (int)floor(fragmentPosition / mGroupSize);
        int groupPosEnd = groupPosStart + mGroupSize - 1;

        ReadPositionGroup group = new ReadPositionGroup(groupPosStart, groupPosEnd);
        mPositionGroups.add(groupIndex, group);
        return group;
    }

    private class ReadPositionGroup
    {
        public final int PositionStart;
        public final int PositionEnd;

        public final Map<FragmentCoords,List<SAMRecord>> FragCoordsMap;

        public ReadPositionGroup(final int positionStart, final int positionEnd)
        {
            PositionStart = positionStart;
            PositionEnd = positionEnd;
            FragCoordsMap = Maps.newHashMap();
        }

        public void addRead(final FragmentCoords fragCoords, final SAMRecord read)
        {
            List<SAMRecord> reads = FragCoordsMap.get(fragCoords.Key);

            if(reads == null)
            {
                reads = Lists.newArrayList();
                FragCoordsMap.put(fragCoords, reads);
            }

            reads.add(read);
        }

        public int readCount() { return FragCoordsMap.values().stream().mapToInt(x -> x.size()).sum(); }

        public String toString()
        {
            return format("range(%d-%d) frags(%d) reads(%d)", PositionStart, PositionEnd, FragCoordsMap.size(), readCount());
        }
    }

    public String toString()
    {
        int fragCount = mPositionGroups.stream().mapToInt(x -> x.FragCoordsMap.size()).sum();
        int readCount = mPositionGroups.stream().mapToInt(x -> x.readCount()).sum();

        return format("minPos(%s:d) groups(%d) frags(%d) reads(%d)",
                mCurrentChromosome, mCurrentReadMinPosition, mPositionGroups.size(), fragCount, readCount);
    }
}
