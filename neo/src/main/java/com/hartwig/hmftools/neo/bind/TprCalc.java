package com.hartwig.hmftools.neo.bind;

import java.util.List;

import com.google.common.collect.Lists;

public class TprCalc
{
    private final int[] mCounts;
    private int mTotal;

    private static final List<Double> RANK_BUCKETS = Lists.newArrayList(0.00025, 0.0005, 0.001, 0.002, 0.005);

    public TprCalc()
    {
        mCounts = new int[RANK_BUCKETS.size()];
        mTotal = 0;
    }

    public int entryCount() { return mTotal; }

    public void addRank(double rank)
    {
        ++mTotal;

        for(int i = 0; i < RANK_BUCKETS.size(); ++i)
        {
            if(rank <= RANK_BUCKETS.get(i))
                ++mCounts[i];
        }
    }

    public double calc()
    {
        double meanTotal = 0;

        for(int i = 0; i < RANK_BUCKETS.size(); ++i)
        {
            double mean = mCounts[i] / (double)mTotal;
            meanTotal += mean;
        }

        return meanTotal / RANK_BUCKETS.size();
    }
}
