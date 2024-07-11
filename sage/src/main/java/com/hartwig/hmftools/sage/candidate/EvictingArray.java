package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.function.Consumer;
import java.util.function.Function;

import com.google.common.annotations.VisibleForTesting;

public class EvictingArray
{
    // a ring buffer to store alts and ref/alt depth
    // capacity is designed to cover a set of reads covering at least a few multiples of the maximum read length
    private final RefContext[] mElements;
    private final Consumer<RefContext> mEvictionHandler;
    private int mMinPosition;
    private int mMinPositionIndex;

    private final int mReadLengthBuffer;
    private final int mCapacity;

    public static final int MAX_EXPECTED_DEL = 80;

    public EvictingArray(int maxReadLength, final Consumer<RefContext> evictionHandler)
    {
        mEvictionHandler = evictionHandler;
        mReadLengthBuffer = maxReadLength + MAX_EXPECTED_DEL;
        mCapacity = maxReadLength * 2;
        mElements = new RefContext[mCapacity];
        mMinPosition = 0;
        mMinPositionIndex = 0;
    }

    public RefContext getOrCreateRefContext(int position, final Function<Integer,RefContext> supplier)
    {
        if(!isValidPosition(position))
            return null;

        int minPositionStart = mMinPosition;
        checkFlush(position);

        int distanceFromMinPosition = position - mMinPosition;
        int index = calcIndex(distanceFromMinPosition);

        if(index < 0 || index >= mCapacity)
        {
            SG_LOGGER.error("invalid evicting array index({}) position({}) minPosition({} -> {})",
                    index, position, minPositionStart, mMinPosition);
            System.exit(1);
        }

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

    private boolean isValidPosition(int position)
    {
        if(mMinPosition > 0 && position < mMinPosition)
        {
            SG_LOGGER.warn("ignoring read with position({}) before prior position({})", position, mMinPosition);
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

            // move min position to at most the read length behind the current position
            int newMinPosition = position - mReadLengthBuffer;
            flushCount = min(newMinPosition - mMinPosition, mCapacity);
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
        mMinPosition = max(1, position - mReadLengthBuffer);
    }

    @VisibleForTesting
    public int minPosition() { return mMinPosition; }

    @VisibleForTesting
    public int readLengthBuffer() { return mReadLengthBuffer; }

    @VisibleForTesting
    public int capacity() { return mCapacity; }

    @VisibleForTesting
    public int itemCount()
    {
        int items = 0;

        for(int i = 0; i < mElements.length; ++i)
        {
            if(mElements[i] != null)
                ++items;
        }

        return items;
    }

    public String toString()
    {
        return format("capacity(%d) buffer(%d) minPos(%d) itemCount(%d)",
                mCapacity, mReadLengthBuffer, mMinPosition, itemCount());
    }
}
