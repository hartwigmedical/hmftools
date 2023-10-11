package com.hartwig.hmftools.bamtools.metrics;

public class FlagQCStats
{
    private int mPassed;
    private int mFailed;

    public FlagQCStats()
    {
        mPassed = 0;
        mFailed = 0;
    }

    public void record(boolean passes)
    {
        if(passes)
        {
            mPassed++;
        }
        else
        {
            mFailed++;
        }
    }

    public void merge(final FlagQCStats other)
    {
        mPassed += other.mPassed;
        mFailed += other.mFailed;
    }

    public int getPassed()
    {
        return mPassed;
    }

    public int getFailed()
    {
        return mFailed;
    }

    @Override
    public String toString()
    {
        return String.format("%d + %d", mPassed, mFailed);
    }
}
