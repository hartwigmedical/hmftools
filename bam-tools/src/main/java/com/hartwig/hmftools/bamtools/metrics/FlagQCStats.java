package com.hartwig.hmftools.bamtools.metrics;

public class FlagQCStats
{
    private long mPassed;
    private long mFailed;

    public FlagQCStats()
    {
        mPassed = 0;
        mFailed = 0;
    }

    public void record(boolean passes)
    {
        if(passes)
            mPassed++;
        else
            mFailed++;
    }

    public void decrement(boolean passes)
    {
        if(passes)
            --mPassed;
        else
            --mFailed;
    }

    public void merge(final FlagQCStats other)
    {
        mPassed += other.mPassed;
        mFailed += other.mFailed;
    }

    public long getPassed()
    {
        return mPassed;
    }
    public long getFailed()
    {
        return mFailed;
    }

    @Override
    public String toString()
    {
        return String.format("%d + %d", mPassed, mFailed);
    }

    public static String flagStatsPercentages(final FlagQCStats numerator, final FlagQCStats denominator)
    {
        String passedStr = (denominator.getPassed() == 0) ?
                "N/A" : String.format("%.2f%%", 100.0f * numerator.getPassed() / denominator.getPassed());

        String failedStr = (denominator.getFailed() == 0) ?
                "N/A" : String.format("%.2f%%", 100.0f * numerator.getFailed() / denominator.getFailed());

        return String.format("(%s : %s)", passedStr, failedStr);
    }
}
