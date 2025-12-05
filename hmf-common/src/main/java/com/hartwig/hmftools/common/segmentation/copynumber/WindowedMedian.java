package com.hartwig.hmftools.common.segmentation.copynumber;

import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * Efficiently computes running medians using two heaps.
 */
class WindowedMedian
{
    // maxHeap contains the smaller half of elements (including the median if windowSize is odd)
    // minHeap contains the larger half of elements
    private final PriorityQueue<Double> maxHeap;
    private final PriorityQueue<Double> minHeap;
    private final double[] result;

    WindowedMedian(double[] data, int windowSize)
    {
        if(windowSize % 2 != 1)
        {
            throw new IllegalArgumentException("Window size must be odd");
        }

        this.maxHeap = new PriorityQueue<>(Comparator.reverseOrder());
        this.minHeap = new PriorityQueue<>();
        this.result = new double[data.length];

        int halfWindow = windowSize / 2;

        // Initialize heaps with the first window
        for(int i = 0; i < windowSize && i < data.length; i++)
        {
            addToHeaps(data[i]);
        }

        // Calculate median for the first interior position
        if(halfWindow < data.length)
        {
            result[halfWindow] = getMedian();
        }

        // Process the rest of the array
        for(int i = windowSize; i < data.length; i++)
        {
            // Remove the element that's no longer in the window
            removeFromHeaps(data[i - windowSize]);

            // Add the new element
            addToHeaps(data[i]);

            // Calculate median for the current position
            result[i - halfWindow] = getMedian();
        }

        // Fill in the edges (where we don't have a full window)
        for(int i = 0; i < halfWindow && i < data.length; i++)
        {
            result[i] = data[i];
        }
        if(data.length - (data.length - halfWindow) >= 0)
        {
            System.arraycopy(data, data.length - halfWindow, result, data.length - halfWindow, data.length - (data.length - halfWindow));
        }
    }

    double[] getMedians()
    {
        return result;
    }

    private void addToHeaps(double value)
    {
        if(maxHeap.isEmpty() || value <= maxHeap.peek())
        {
            maxHeap.add(value);
        }
        else
        {
            minHeap.add(value);
        }
        balanceHeaps();
    }

    private void removeFromHeaps(double value)
    {
        if(value <= maxHeap.peek())
        {
            maxHeap.remove(value);
        }
        else
        {
            minHeap.remove(value);
        }
        balanceHeaps();
    }

    private void balanceHeaps()
    {
        if(maxHeap.size() > minHeap.size() + 1)
        {
            minHeap.add(maxHeap.poll());
        }
        else if(maxHeap.size() < minHeap.size())
        {
            maxHeap.add(minHeap.poll());
        }
    }

    // The number of elements in the heaps is odd so the median is at the top of the max heap.
    private double getMedian()
    {
        return maxHeap.peek();
    }
}