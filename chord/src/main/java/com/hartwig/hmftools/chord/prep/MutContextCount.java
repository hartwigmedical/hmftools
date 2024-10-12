package com.hartwig.hmftools.chord.prep;

import java.util.Objects;

public class MutContextCount
{
    public final String mName;
    public final int mCount;

    public MutContextCount(String contextName, int count)
    {
        mName = contextName;
        mCount = count;
    }

    @Override
    public String toString()
    {
        return mName + "\t" + mCount;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mName, mCount);
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

        final MutContextCount other = (MutContextCount) o;
        return mCount == other.mCount && mName.equals(other.mName);
    }
}
