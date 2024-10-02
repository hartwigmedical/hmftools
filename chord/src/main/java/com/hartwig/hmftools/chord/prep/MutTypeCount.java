package com.hartwig.hmftools.chord.prep;

public class MutTypeCount
{
    String mType;
    int mCount;

    public MutTypeCount(String mutType, int count)
    {
        mType = mutType;
        mCount = count;
    }

    @Override
    public String toString()
    {
        return mType + "\t" + mCount;
    }
}
