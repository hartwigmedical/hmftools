package com.hartwig.hmftools.sieve.bedstats;

public class Stats
{
    public static final String TSV_HEADER = "MaxDepth\tAvgDepth";

    private final int mBaseLength;
    private int mMaxDepth;
    private long mTotalVolume;

    public Stats(final int[] baseDepth)
    {
        mBaseLength = baseDepth.length;
        mMaxDepth = 0;
        mTotalVolume = 0;
        for(int i = 0; i < baseDepth.length; i++)
        {
            mMaxDepth = Math.max(mMaxDepth, baseDepth[i]);
            mTotalVolume += baseDepth[i];
        }
    }

    public String getTSVFragment()
    {
        final float avgDepth = 1.0f * mTotalVolume / mBaseLength;
        return String.valueOf(mMaxDepth) + '\t' + avgDepth;
    }
}
