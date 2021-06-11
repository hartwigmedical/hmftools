package com.hartwig.hmftools.sage.quality;

import java.util.concurrent.atomic.AtomicInteger;

import org.jetbrains.annotations.NotNull;

class QualityCounter implements QualityCounterKey, Comparable<QualityCounter>
{
    private final QualityCounterKey mKey;
    private final AtomicInteger mCount = new AtomicInteger();

    public QualityCounter(final QualityCounterKey key)
    {
        mKey = key;
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

    public int position()
    {
        return mKey.position();
    }

    public byte ref()
    {
        return mKey.ref();
    }

    public byte alt()
    {
        return mKey.alt();
    }

    public byte qual()
    {
        return mKey.qual();
    }

    @Override
    public byte[] trinucleotideContext()
    {
        return mKey.trinucleotideContext();
    }

    @NotNull
    public QualityCounterKey key()
    {
        return mKey;
    }

    @Override
    public int compareTo(@NotNull final QualityCounter o2)
    {
        int countCompare = Integer.compare(o2.count(), count());
        if(countCompare != 0)
        {
            return countCompare;
        }

        int refCompare = Byte.compare(ref(), o2.ref());
        if(refCompare != 0)
        {
            return refCompare;
        }

        int altCompare = Byte.compare(alt(), o2.alt());
        if(altCompare != 0)
        {
            return altCompare;
        }

        if(trinucleotideContext().length < 3 || o2.trinucleotideContext().length < 3)
        {
            return 0;
        }

        int triOne = Byte.compare(trinucleotideContext()[0], o2.trinucleotideContext()[0]);
        if(triOne != 0)
        {
            return triOne;
        }

        int triTwo = Byte.compare(trinucleotideContext()[1], o2.trinucleotideContext()[1]);
        if(triTwo != 0)
        {
            return triTwo;
        }

        int triThree = Byte.compare(trinucleotideContext()[2], o2.trinucleotideContext()[2]);
        if(triThree != 0)
        {
            return triThree;
        }

        return 0;
    }
}
