package com.hartwig.hmftools.sage.bqr;

public class BqrKeyCounter
{
    public final BqrKey Key;

    private int mCount;

    public BqrKeyCounter(final BqrKey key)
    {
        Key = key;
        mCount = 0;
    }

    public int count()
    {
        return mCount;
    }
    public void increment(int increment) { mCount += increment;}

    public byte ref() { return Key.Ref; }
    public byte alt()
    {
        return Key.Alt;
    }
}
