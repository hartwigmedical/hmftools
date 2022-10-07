package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.function.Consumer;
import java.util.function.Function;

import javax.annotation.Nullable;

import com.hartwig.hmftools.sage.candidate.RefContext;

public class EvictingArray
{
    // a ring buffer to store alts and ref/alt depth
    // capacity is designed to cover a set of reads covering at least a few multiples of the maximum read length
    private final RefContext[] mElements;
    private final int[] mDepth;
    private final int[] mDepthLimit;
    private final Consumer<RefContext> mEvictionHandler;
    private int mMinPosition;
    private int mMinPositionIndex;

    private final int mCapacity;

    public static final int MIN_CAPACITY = 256;

    public EvictingArray(int capacity, Consumer<RefContext> evictionHandler)
    {
        mEvictionHandler = evictionHandler;
        mCapacity = capacity;
        mElements = new RefContext[mCapacity];
        mDepth = new int[mCapacity];
        mDepthLimit = new int[mCapacity];
        mMinPosition = 0;
        mMinPositionIndex = 0;
    }

    public int minPosition() { return mMinPosition; }
    public int capacity() { return mCapacity; }

    public Integer getDepth(int position)
    {
        int distanceFromMinPosition = position - mMinPosition;
        if(distanceFromMinPosition < 0 || distanceFromMinPosition >= mCapacity)
            return null;

        int index = calcIndex(distanceFromMinPosition);
        return mDepth[index];
    }

    public @Nullable Boolean exceedsDepthLimit(int position)
    {
        int distanceFromMinPosition = position - mMinPosition;
        if(distanceFromMinPosition < 0 || distanceFromMinPosition >= mCapacity)
            return null;

        int index = calcIndex(distanceFromMinPosition);
        if(mDepthLimit[index] == 0)
            return null;

        return mDepth[index] >= mDepthLimit[index];
    }

    public void registerDepthLimit(int position, int limit)
    {
        if(!isValidPosition(position, "registerDepthLimit"))
            return;

        // find position, flush if required, add depth
        checkFlush(position);

        int distanceFromMinPosition = position - mMinPosition;
        int index = calcIndex(distanceFromMinPosition);
        mDepthLimit[index] = limit;
    }

    public void registerDepth(int position)
    {
        if(!isValidPosition(position, "registerDepth"))
            return;

        checkFlush(position);

        int distanceFromMinPosition = position - mMinPosition;
        int index = calcIndex(distanceFromMinPosition);

        ++mDepth[index];
    }

    public RefContext getOrCreateRefContext(int position, final Function<Integer,RefContext> supplier)
    {
        if(!isValidPosition(position, "getOrCreateRefContext"))
            return null;

        checkFlush(position);

        int distanceFromMinPosition = position - mMinPosition;
        int index = calcIndex(distanceFromMinPosition);

        RefContext element = mElements[index];
        if(element == null)
        {
            element = supplier.apply(position);
            mElements[index] = element;
        }

        return element;
    }

    public void evictAll()
    {
        checkFlush(-1);
    }

    private int calcIndex(int distanceFromMinPosition)
    {
        // capacity = 10, min position = 1, min index = 0, position of 10 is index 9
        // capacity = 10, min position = 2, min index = 1, position of 10 is index 9, position of 11 is 0
        if(mMinPositionIndex + distanceFromMinPosition < mCapacity)
            return mMinPositionIndex + distanceFromMinPosition;

        // index is from start of ring buffer
        return distanceFromMinPosition + mMinPositionIndex - mCapacity;
    }

    private boolean isValidPosition(int position, final String caller)
    {
        if(mMinPosition > 0 && position < mMinPosition)
        {
            SG_LOGGER.warn("{}: ignoring read with position({}) before prior position({})", caller, position, mMinPosition);
            return false;
        }

        return true;
    }

    private void checkFlush(int position)
    {
        int flushCount = 0;

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
            RefContext element = mElements[mMinPositionIndex];

            // clear and process each element and depth
            if(element != null)
            {
                mEvictionHandler.accept(element);
                mElements[mMinPositionIndex] = null;
            }

            mDepth[mMinPositionIndex] = 0;
            mDepthLimit[mMinPositionIndex] = 0;

            mMinPosition++;

            if(mMinPositionIndex + 1 >= mElements.length)
                mMinPositionIndex = 0;
            else
                ++mMinPositionIndex;
        }

        if(flushCount >= mCapacity)
            resetMinPosition(position);
    }

    private void resetMinPosition(int position)
    {
        mMinPositionIndex = 0;
        mMinPosition = max(1, position - (int)round(mCapacity * 0.5));
    }
}
