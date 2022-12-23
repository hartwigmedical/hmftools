package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import htsjdk.samtools.SAMRecord;

public class ReadPositionArray
{
    // a ring buffer to store reads at each read starting position
    private final String mChromosome;
    private final PositionFragments[] mForwardPositionGroups;
    private final Map<Integer, PositionFragments> mReversePositionGroups;
    private final Map<String,Fragment> mFragments;
    private final Consumer<PositionFragments> mReadGroupHandler;
    private int mMinPosition;
    private int mMinPositionIndex;

    private final int mCapacity;

    public ReadPositionArray(final String chromosome, int capacity, final Consumer<PositionFragments> evictionHandler)
    {
        mChromosome = chromosome;
        mReadGroupHandler = evictionHandler;
        mCapacity = capacity;
        mForwardPositionGroups = new PositionFragments[mCapacity];
        mReversePositionGroups = Maps.newHashMap();
        mFragments = Maps.newHashMap();
        mMinPosition = 0;
        mMinPositionIndex = 0;
    }

    public boolean processRead(final SAMRecord read)
    {
        int readPosStart = read.getAlignmentStart();

        if(mMinPosition == 0)
            mMinPosition = readPosStart;

        if(!isValidPosition(readPosStart))
            return false;

        checkFlush(readPosStart);

        boolean readStored = handleRead(read);

        return readStored;
    }

    private boolean handleRead(final SAMRecord read)
    {
        /* scenarios:
            - unpaired or mate is unmapped
                - store as potential duplicate
                - supplementaries - store if within range, to avoid processing later
            - read's mate is on another chromosome - store since can still evaluate if duplicate (from both chromosomes)
            - lower pos with mate higher on same chromosome - store regardless
            - higher pos with mate lower on same chromosome
                - mate within current range - look up and add
                - mate already flushed, return unhandled
            - supplementary
                - with mate within current range - look up and add
                - otherwise return unhandled
        */

        // for supplementaries, add the read if the non-supplmentaries have already been stored
        if(read.getSupplementaryAlignmentFlag())
            return storeAdditionalRead(read);

        if(!read.getReadPairedFlag() || read.getMateUnmappedFlag())
        {
            storeInitialRead(read);
            return true;
        }

        // mate is elsewhere so store this primary read
        if(!read.getMateReferenceName().equals(mChromosome))
        {
            storeInitialRead(read);
            return true;
        }

        Fragment readGroup = mFragments.get(read.getReadName());

        if(readGroup != null)
        {
            readGroup.addRead(read);
            return true;
        }

        storeInitialRead(read);
        return true;
    }

    private void storeInitialRead(final SAMRecord read)
    {
        Fragment fragment = new Fragment(read);
        int fragmentPosition = fragment.initialPosition();

        mFragments.put(read.getReadName(), fragment);

        if(orientation(read) == POS_ORIENT)
        {
            int index = calcIndex(fragmentPosition);

            PositionFragments element = mForwardPositionGroups[index];
            if(element == null)
            {
                element = new PositionFragments(fragment, fragmentPosition);
                mForwardPositionGroups[index] = element;
            }
            else
            {
                element.Fragments.add(fragment);
            }
        }
        else
        {
            // store in reverse strand map
            PositionFragments element = mReversePositionGroups.get(fragmentPosition);
            if(element == null)
            {
                element = new PositionFragments(fragment, fragmentPosition);
                mReversePositionGroups.put(fragmentPosition, element);
            }

            element.Fragments.add(fragment);
        }
    }

    private boolean storeAdditionalRead(final SAMRecord read)
    {
        Fragment fragment = mFragments.get(read.getReadName());

        if(fragment != null)
        {
            fragment.addRead(read);
            return true;
        }

        return false;
    }

    public void evictAll()
    {
        checkFlush(-1);
    }

    private int calcIndex(int position)
    {
        // capacity = 10, min position = 1, min index = 0, position of 10 is index 9
        // capacity = 10, min position = 2, min index = 1, position of 10 is index 9, position of 11 is 0
        int distanceFromMinPosition = position - mMinPosition;

        if(mMinPositionIndex + distanceFromMinPosition < mCapacity)
            return mMinPositionIndex + distanceFromMinPosition;

        // index is from start of ring buffer
        return distanceFromMinPosition + mMinPositionIndex - mCapacity;
    }

    private boolean isValidPosition(int position)
    {
        if(mMinPosition > 0 && position < mMinPosition)
        {
            BM_LOGGER.warn("ignoring read with position({}) before prior position({})", position, mMinPosition);
            return false;
        }

        return true;
    }

    private void checkFlush(int position)
    {
        int flushCount = 0;

        Set<String> flushedReadIds = null;

        if(position > 0)
        {
            if(mMinPosition == 0)
            {
                resetMinPosition(position);
                return;
            }

            int distanceFromMinPosition = position - mMinPosition;

            if(distanceFromMinPosition < mCapacity)
                return;

            flushCount = position - mMinPosition - mCapacity + 1;
        }
        else
        {
            flushCount = mCapacity;
        }

        // only iterate at most once through the array
        for(int i = 0; i < min(flushCount, mCapacity); i++)
        {
            PositionFragments element = mForwardPositionGroups[mMinPositionIndex];

            // clear and process each element and depth
            if(element != null)
            {
                if(flushedReadIds == null)
                    flushedReadIds = Sets.newHashSet();

                mReadGroupHandler.accept(element);

                for(Fragment fragment : element.Fragments)
                    flushedReadIds.add(fragment.id());

                mForwardPositionGroups[mMinPositionIndex].Fragments.clear();
                mForwardPositionGroups[mMinPositionIndex] = null;
            }

            mMinPosition++;

            if(mMinPositionIndex + 1 >= mForwardPositionGroups.length)
                mMinPositionIndex = 0;
            else
                ++mMinPositionIndex;
        }

        if(flushCount >= mCapacity)
            resetMinPosition(position);

        if(!flushedReadIds.isEmpty())
        {
            Set<Integer> flushedPositions = Sets.newHashSet();
            for(Map.Entry<Integer, PositionFragments> entry : mReversePositionGroups.entrySet())
            {
                if(position < 0 || entry.getKey() < position)
                {
                    flushedPositions.add(entry.getKey());

                    for(Fragment fragment : entry.getValue().Fragments)
                        flushedReadIds.add(fragment.id());

                    mReadGroupHandler.accept(entry.getValue());
                }
            }

            flushedPositions.forEach(x -> mReversePositionGroups.remove(x));
        }
    }

    private void resetMinPosition(int position)
    {
        mMinPositionIndex = 0;
        mMinPosition = max(1, position - (int)round(mCapacity * 0.5));
    }
}
