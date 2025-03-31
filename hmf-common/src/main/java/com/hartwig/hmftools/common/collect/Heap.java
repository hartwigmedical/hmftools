package com.hartwig.hmftools.common.collect;

import java.util.Comparator;
import java.util.PriorityQueue;

public class Heap<T>
{
    private class Element implements Comparable<Element>
    {
        public final T Element;
        public final long Count;

        public Element(final T element)
        {
            Element = element;
            Count = mCounter;
            mCounter++;
        }

        @Override
        public int compareTo(final Heap<T>.Element o)
        {
            int elementComparison = mComparator.compare(Element, o.Element);
            if(elementComparison != 0)
                return elementComparison;

            if(Count < o.Count)
                return -1;

            if(Count > o.Count)
                return 1;

            return 0;
        }
    }

    private final Comparator<T> mComparator;
    private final PriorityQueue<Element> mQueue;
    private long mCounter;

    public Heap(final Comparator<T> comparator)
    {
        mComparator = comparator;
        mQueue = new PriorityQueue<>();
        mCounter = 0;
    }

    public int size()
    {
        return mQueue.size();
    }

    public boolean isEmpty()
    {
        return mQueue.isEmpty();
    }

    public void clear()
    {
        mQueue.clear();
    }

    public void add(final T element)
    {
        mQueue.add(new Element(element));
    }

    public void addAll(final Iterable<T> elements)
    {
        for(T el : elements)
            mQueue.add(new Element(el));
    }

    public T pop()
    {
        Element el = mQueue.poll();
        return el == null ? null : el.Element;
    }

    public T peek()
    {
        Element el = mQueue.peek();
        return el == null ? null : el.Element;
    }
}