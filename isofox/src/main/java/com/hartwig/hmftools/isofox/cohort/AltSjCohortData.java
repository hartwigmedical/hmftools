package com.hartwig.hmftools.isofox.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

public class AltSjCohortData
{
    public final AltSpliceJunction AltSJ;

    // cohort data
    private final List<String> mSampleIdsCohortA;
    private final List<String> mSampleIdsCohortB;

    private int mTotalFragmentCountCohortA;
    private int mTotalFragmentCountCohortB;

    private int mMaxFragmentCountCohortA;
    private int mMaxFragmentCountCohortB;
    private final int[] mPositionCounts; // counts at the start and end

    public AltSjCohortData(final AltSpliceJunction altSJ)
    {
        AltSJ = altSJ;
        mSampleIdsCohortA = Lists.newArrayList();
        mSampleIdsCohortB = Lists.newArrayList();

        mTotalFragmentCountCohortA = 0;
        mTotalFragmentCountCohortB = 0;

        mMaxFragmentCountCohortA = 0;
        mMaxFragmentCountCohortB = 0;
        mPositionCounts = new int[SE_PAIR];
    }

    public void addSampleAndCount(final String sampleId, int fragCount, boolean isCohortA)
    {
        if(isCohortA)
        {
            mSampleIdsCohortA.add(sampleId);
            mTotalFragmentCountCohortA += fragCount;
            mMaxFragmentCountCohortA = max(mMaxFragmentCountCohortA, fragCount);
        }
        else
        {
            mSampleIdsCohortB.add(sampleId);
            mTotalFragmentCountCohortB += fragCount;
            mMaxFragmentCountCohortB = max(mMaxFragmentCountCohortB, fragCount);
        }
    }

    public int totalSamples() { return mSampleIdsCohortA.size() + mSampleIdsCohortB.size(); }

    public final List<String> getSampleIds(boolean isCohortA) { return isCohortA ? mSampleIdsCohortA : mSampleIdsCohortB; }

    public int getMaxFragmentCount(boolean isCohortA) { return isCohortA ? mMaxFragmentCountCohortA : mMaxFragmentCountCohortB; }

    public double getAvgFragmentCount(boolean isCohortA)
    {
        return isCohortA ?
                (mSampleIdsCohortA.size() > 0 ? mTotalFragmentCountCohortA / (double)mSampleIdsCohortA.size() : 0)
                : (mSampleIdsCohortB.size() > 0 ? mTotalFragmentCountCohortB / (double)mSampleIdsCohortB.size() : 0);
    }

    public int getPositionCount(int seIndex) { return mPositionCounts[seIndex]; }
    public void addPositionCount(int seIndex, int count) { mPositionCounts[seIndex] += count; }



}
