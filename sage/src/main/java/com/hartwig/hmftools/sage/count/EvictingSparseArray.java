package com.hartwig.hmftools.sage.count;

import java.util.function.Consumer;
import java.util.function.Function;

public class EvictingSparseArray<T> {

    private final Object[] elements;
    private final Consumer<T> evictionHandler;
    private int minPosition = 0;
    private int minPositionIndex = 0;

    private final int maxSize;

    public EvictingSparseArray(int maxSize, Consumer<T> evictionHandler) {
        this.evictionHandler = evictionHandler;
        this.elements = new Object[maxSize];
        this.maxSize = maxSize;
    }

    public T computeIfAbsent(long position, Function<Long, T> supplier) {
        if (minPosition == 0) {
            minPosition = (int) position - maxSize + 1;
        }

        int distanceFromMinPosition = (int) position - minPosition;
        if (distanceFromMinPosition < 0) {
            throw new IllegalArgumentException("Cannot add position: " + position + " before min position: " + minPosition);
        }

        if (distanceFromMinPosition >= maxSize) {
            flush((int) position - minPosition - maxSize + 1);
        }

        distanceFromMinPosition = (int) position - minPosition;
        int index = (minPositionIndex + distanceFromMinPosition) & (elements.length - 1);
        T element = (T) elements[index];
        if (element == null) {
            element = supplier.apply(position);
            elements[index] = element;
        }

        return element;
    }


    int minPosition() {
        return minPosition;
    }

    public void evictAll() {
        flush(maxSize);
    }

    private void flush(int count) {
        for (int i = 0; i < count; i++) {
            T element = (T) elements[minPositionIndex];
            if (element != null) {
                evictionHandler.accept(element);
                elements[minPositionIndex] = null;
            }
            minPosition++;
            minPositionIndex = (minPositionIndex + 1) & (elements.length - 1);
        }
    }

}
