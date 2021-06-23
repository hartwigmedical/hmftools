package com.hartwig.hmftools.sage.quality;

import java.util.concurrent.atomic.AtomicInteger;

class QualityCounter implements Comparable<QualityCounter>
{
    public final QualityCounterKey Key;

    private final AtomicInteger mCount = new AtomicInteger();

    public QualityCounter(final QualityCounterKey key)
    {
        Key = key;
    }

    public int count()
    {
        return mCount.get();
    }

    public void increment()
    {
        increment(1);
    }
    public void increment(int increment)
    {
        mCount.addAndGet(increment);
    }

    public int position() { return Key.Position; }

    public byte ref() { return Key.Ref; }
    public byte alt()
    {
        return Key.Alt;
    }
    public byte qual()
    {
        return Key.Quality;
    }

    public byte[] trinucleotideContext()
    {
        return Key.TrinucleotideContext;
    }

    @Override
    public int compareTo(final QualityCounter other)
    {
        int countCompare = Integer.compare(other.count(), count());
        if(countCompare != 0)
            return countCompare;

        return Key.compareTo(other.Key, false);
    }
}
