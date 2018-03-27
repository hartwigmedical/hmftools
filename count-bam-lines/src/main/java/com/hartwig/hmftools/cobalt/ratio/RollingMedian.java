package com.hartwig.hmftools.cobalt.ratio;

import java.util.Comparator;
import java.util.PriorityQueue;

class RollingMedian {

    private final PriorityQueue<Double> minheap = new PriorityQueue<>(new MinHeapComparator());
    private final PriorityQueue<Double> maxheap = new PriorityQueue<>(new MaxHeapComparator());

    void add(double n) {
        if (isEmpty()) {
            minheap.add(n);
        } else {
            if (Double.compare(n, median()) <= 0) {
                maxheap.add(n);
            } else {
                minheap.add(n);
            }
        }
        fixChaos();
    }

    void remove(double n) {
        if (!isEmpty()) {
            if (Double.compare(n, median()) <= 0) {
                maxheap.remove(n);
            } else {
                minheap.remove(n);
            }
        }
        fixChaos();
    }

    int size() {
        return maxheap.size() + minheap.size();
    }

    private boolean isEmpty() {
        return size() == 0;
    }

    private void fixChaos() {
        //if sizes of heaps differ by 2, then it's a chaos, since median must be the middle element
        if (Math.abs(maxheap.size() - minheap.size()) > 1) {
            //check which one is the culprit and take action by kicking out the root from culprit into victim
            if (maxheap.size() > minheap.size()) {
                minheap.add(maxheap.poll());
            } else {
                maxheap.add(minheap.poll());
            }
        }
    }

    double median() {
        if (isEmpty()) {
            return 0;
        }
        if (maxheap.size() == minheap.size()) {
            return (maxheap.peek() + minheap.peek()) / 2;
        } else if (maxheap.size() > minheap.size()) {
            return maxheap.peek();
        } else {
            return minheap.peek();
        }
    }

    private static class MinHeapComparator implements Comparator<Double> {
        @Override
        public int compare(Double i, Double j) {
            return Double.compare(i, j);
        }
    }

    private static class MaxHeapComparator implements Comparator<Double> {
        // opposite to minHeapComparator, invert the return values
        @Override
        public int compare(Double i, Double j) {
            return -1 * Double.compare(i, j);
        }
    }
}