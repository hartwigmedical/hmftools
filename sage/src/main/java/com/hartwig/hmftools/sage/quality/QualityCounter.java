package com.hartwig.hmftools.sage.quality;

import java.util.concurrent.atomic.AtomicInteger;

import org.jetbrains.annotations.NotNull;

class QualityCounter implements QualityCounterKey, Comparable<QualityCounter>
{

    private final QualityCounterKey key;
    private final AtomicInteger count = new AtomicInteger();

    public QualityCounter(final QualityCounterKey key)
    {
        this.key = key;
    }

    public int count()
    {
        return count.get();
    }

    public void increment()
    {
        increment(1);
    }

    public void increment(int increment)
    {
        count.addAndGet(increment);
    }

    public int position()
    {
        return key.position();
    }

    public byte ref()
    {
        return key.ref();
    }

    public byte alt()
    {
        return key.alt();
    }

    public byte qual()
    {
        return key.qual();
    }

    @Override
    public byte[] trinucleotideContext()
    {
        return key.trinucleotideContext();
    }

    @NotNull
    public QualityCounterKey key()
    {
        return key;
    }

    @Override
    public int compareTo(@NotNull final QualityCounter o2)
    {
        int countCompare = Integer.compare(o2.count(), this.count());
        if(countCompare != 0)
        {
            return countCompare;
        }

        int refCompare = Byte.compare(this.ref(), o2.ref());
        if(refCompare != 0)
        {
            return refCompare;
        }

        int altCompare = Byte.compare(this.alt(), o2.alt());
        if(altCompare != 0)
        {
            return altCompare;
        }

        if(this.trinucleotideContext().length < 3 || o2.trinucleotideContext().length < 3)
        {
            return 0;
        }

        int triOne = Byte.compare(this.trinucleotideContext()[0], o2.trinucleotideContext()[0]);
        if(triOne != 0)
        {
            return triOne;
        }

        int triTwo = Byte.compare(this.trinucleotideContext()[1], o2.trinucleotideContext()[1]);
        if(triTwo != 0)
        {
            return triTwo;
        }

        int triThree = Byte.compare(this.trinucleotideContext()[2], o2.trinucleotideContext()[2]);
        if(triThree != 0)
        {
            return triThree;
        }

        return 0;
    }
}
