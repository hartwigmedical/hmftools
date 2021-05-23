package com.hartwig.hmftools.cobalt.ratio;

import java.util.Comparator;
import java.util.PriorityQueue;

public class RollingMedian
{
    private final PriorityQueue<Double> mMinheap;
    private final PriorityQueue<Double> mMaxheap;

    public RollingMedian()
    {
       mMinheap = new PriorityQueue<>(new MinHeapComparator());
       mMaxheap = new PriorityQueue<>(new MaxHeapComparator());
    }

    public void add(double n)
    {
        if(isEmpty())
        {
            mMinheap.add(n);
        }
        else
        {
            if(Double.compare(n, median()) <= 0)
            {
                mMaxheap.add(n);
            }
            else
            {
                mMinheap.add(n);
            }
        }

        fixChaos();
    }

    public void remove(double n)
    {
        if(!isEmpty())
        {
            if(Double.compare(n, median()) <= 0)
            {
                mMaxheap.remove(n);
            }
            else
            {
                mMinheap.remove(n);
            }
        }
        fixChaos();
    }

    public int size()
    {
        return mMaxheap.size() + mMinheap.size();
    }

    private boolean isEmpty()
    {
        return size() == 0;
    }

    private void fixChaos()
    {
        //if sizes of heaps differ by 2, then it's a chaos, since median must be the middle element
        if(Math.abs(mMaxheap.size() - mMinheap.size()) > 1)
        {
            //check which one is the culprit and take action by kicking out the root from culprit into victim
            if(mMaxheap.size() > mMinheap.size())
            {
                mMinheap.add(mMaxheap.poll());
            }
            else
            {
                mMaxheap.add(mMinheap.poll());
            }
        }
    }

    public double median()
    {
        if(isEmpty())
        {
            return 0;
        }
        if(mMaxheap.size() == mMinheap.size())
        {
            return (mMaxheap.peek() + mMinheap.peek()) / 2;
        }
        else if(mMaxheap.size() > mMinheap.size())
        {
            return mMaxheap.peek();
        }
        else
        {
            return mMinheap.peek();
        }
    }

    private static class MinHeapComparator implements Comparator<Double>
    {
        @Override
        public int compare(Double i, Double j)
        {
            return Double.compare(i, j);
        }
    }

    private static class MaxHeapComparator implements Comparator<Double>
    {
        // opposite to minHeapComparator, invert the return values
        @Override
        public int compare(Double i, Double j)
        {
            return -1 * Double.compare(i, j);
        }
    }
}