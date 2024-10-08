package com.hartwig.hmftools.chord.common;

import java.util.Objects;

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

    @Override
    public int hashCode()
    {
        return Objects.hash(mType, mCount);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }

        if(o == null || getClass() != o.getClass())
        {
            return false;
        }

        final MutTypeCount other = (MutTypeCount) o;
        return mCount == other.mCount && mType.equals(other.mType);
    }
}
