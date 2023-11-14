package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

public class StrandBiasData
{
    private int mForwardCount;
    private int mReverseCount;

    public StrandBiasData()
    {
        mForwardCount = 0;
        mReverseCount = 0;
    }

    public void add(boolean isFoward)
    {
        if(isFoward)
            ++mForwardCount;
        else
            ++mReverseCount;
    }

    public int depth() { return mForwardCount + mReverseCount; }

    public double bias()
    {
        double depth = mForwardCount + mReverseCount;
        return depth > 0 ? mForwardCount / depth : 0.5;
    }

    public String toString() { return format("fwd=%d rev=%d total=%d bias=%.3f", mForwardCount, mReverseCount, depth(), bias()); }
}
