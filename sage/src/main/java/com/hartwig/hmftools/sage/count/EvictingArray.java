package com.hartwig.hmftools.sage.count;

import java.util.function.Consumer;
import java.util.function.Function;

import org.jetbrains.annotations.NotNull;

public class EvictingArray<T>
{

    private final Object[] elements;
    private final Consumer<T> evictionHandler;
    private int minPosition = 0;
    private int minPositionIndex = 0;

    private final int capacity;

    public EvictingArray(int minCapacity, Consumer<T> evictionHandler)
    {
        this.evictionHandler = evictionHandler;
        this.capacity = calculateSize(minCapacity);
        this.elements = new Object[this.capacity];
    }

    public T computeIfAbsent(long position, @NotNull final Function<Long, T> supplier)
    {
        if(minPosition == 0)
        {
            minPosition = (int) position - capacity + 1;
        }

        int distanceFromMinPosition = (int) position - minPosition;
        if(distanceFromMinPosition < 0)
        {
            throw new IllegalArgumentException("Cannot add position: " + position + " before min position: " + minPosition);
        }

        if(distanceFromMinPosition >= capacity)
        {
            flush((int) position - minPosition - capacity + 1);
        }

        distanceFromMinPosition = (int) position - minPosition;
        int index = (minPositionIndex + distanceFromMinPosition) & (elements.length - 1);
        T element = (T) elements[index];
        if(element == null)
        {
            element = supplier.apply(position);
            elements[index] = element;
        }

        return element;
    }

    public int capacity()
    {
        return capacity;
    }

    int minPosition()
    {
        return minPosition;
    }

    public void evictAll()
    {
        flush(capacity);
    }

    private void flush(int count)
    {
        for(int i = 0; i < count; i++)
        {
            T element = (T) elements[minPositionIndex];
            if(element != null)
            {
                evictionHandler.accept(element);
                elements[minPositionIndex] = null;
            }
            minPosition++;
            minPositionIndex = (minPositionIndex + 1) & (elements.length - 1);
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
