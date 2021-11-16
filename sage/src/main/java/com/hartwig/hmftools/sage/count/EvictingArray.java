package com.hartwig.hmftools.sage.count;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.function.Consumer;
import java.util.function.Function;

import org.jetbrains.annotations.NotNull;

public class EvictingArray<T>
{
    private final Object[] mElements;
    private final Consumer<T> mEvictionHandler;
    private int mMinPosition = 0;
    private int mMinPositionIndex = 0;

    private final int mCapacity;

    public EvictingArray(int minCapacity, Consumer<T> evictionHandler)
    {
        mEvictionHandler = evictionHandler;
        mCapacity = calculateSize(minCapacity);
        mElements = new Object[mCapacity];
    }

    public T computeIfAbsent(long position, @NotNull final Function<Long, T> supplier)
    {
        if(mMinPosition == 0)
        {
            mMinPosition = (int) position - mCapacity + 1;
        }

        int distanceFromMinPosition = (int) position - mMinPosition;
        if(distanceFromMinPosition < 0)
        {
            SG_LOGGER.warn("ignoring read with position({}) before prior position({})", position, mMinPosition);
            return null;
            // throw new IllegalArgumentException("Cannot add position: " + position + " before min position: " + mMinPosition);
        }

        if(distanceFromMinPosition >= mCapacity)
        {
            flush((int) position - mMinPosition - mCapacity + 1);
        }

        distanceFromMinPosition = (int) position - mMinPosition;
        int index = (mMinPositionIndex + distanceFromMinPosition) & (mElements.length - 1);
        T element = (T) mElements[index];
        if(element == null)
        {
            element = supplier.apply(position);
            mElements[index] = element;
        }

        return element;
    }

    public int capacity()
    {
        return mCapacity;
    }

    int minPosition()
    {
        return mMinPosition;
    }

    public void evictAll()
    {
        flush(mCapacity);
    }

    private void flush(int count)
    {
        for(int i = 0; i < count; i++)
        {
            T element = (T) mElements[mMinPositionIndex];
            if(element != null)
            {
                mEvictionHandler.accept(element);
                mElements[mMinPositionIndex] = null;
            }
            mMinPosition++;
            mMinPositionIndex = (mMinPositionIndex + 1) & (mElements.length - 1);
        }
    }

    private static final int MIN_INITIAL_CAPACITY = 8;

    private static int calculateSize(int numElements)
    {
        int initialCapacity = MIN_INITIAL_CAPACITY;
        // Find the best power of two to hold elements.
        // Tests "<=" because arrays aren't kept full.
        if(numElements >= initialCapacity)
        {
            initialCapacity = numElements - 1;
            initialCapacity |= (initialCapacity >>> 1);
            initialCapacity |= (initialCapacity >>> 2);
            initialCapacity |= (initialCapacity >>> 4);
            initialCapacity |= (initialCapacity >>> 8);
            initialCapacity |= (initialCapacity >>> 16);
            initialCapacity++;

            if(initialCapacity < 0)   // Too many elements, must back off
            {
                initialCapacity >>>= 1;// Good luck allocating 2 ^ 30 elements
            }
        }
        return initialCapacity;
    }
}
