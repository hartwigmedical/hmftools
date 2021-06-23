package com.hartwig.hmftools.sage.quality;

class QualityCounter
{
    public final BaseQualityKey Key;

    private int mCount;

    public QualityCounter(final BaseQualityKey key)
    {
        Key = key;
        mCount = 0;
    }

    public int count()
    {
        return mCount;
    }

    public void increment() { ++mCount; }
    public void increment(int increment) { mCount += increment;}

    public byte ref() { return Key.Ref; }
    public byte alt()
    {
        return Key.Alt;
    }

    public byte[] trinucleotideContext()
    {
        return Key.TrinucleotideContext;
    }

}
