package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

public class AltSpliceJuncCohortData
{
    public final AltSpliceJunction AltSJ;

    // cohort data
    private final List<String> mSampleIds;
    private final List<String> mSampleIdsCohortA;
    private final List<String> mSampleIdsCohortB;

    private int mTotalFragmentCount;
    private int mTotalFragmentCountCohortA;
    private int mTotalFragmentCountCohortB;

    private int mMaxFragmentCount;
    private int mMaxFragmentCountCohortA;
    private int mMaxFragmentCountCohortB;
    private final int[] mPositionCounts; // counts at the start and end

    public AltSpliceJuncCohortData(final AltSpliceJunction altSJ)
    {
        AltSJ = altSJ;
        mSampleIds = Lists.newArrayList();
        mSampleIdsCohortA = Lists.newArrayList();
        mSampleIdsCohortB = Lists.newArrayList();

        mTotalFragmentCount = 0;
        mTotalFragmentCountCohortA = 0;
        mTotalFragmentCountCohortB = 0;

        mMaxFragmentCount = 0;
        mMaxFragmentCountCohortA = 0;
        mMaxFragmentCountCohortB = 0;
        mPositionCounts = new int[SE_PAIR];
    }

    public void addSampleAndCount(final String sampleId, int fragCount)
    {
        mSampleIds.add(sampleId);
        mTotalFragmentCount += fragCount;
        mMaxFragmentCount = max(mMaxFragmentCount, fragCount);
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

    public final List<String> getSampleIds() { return mSampleIds; }
    public final List<String> getCohortSampleIds(boolean isCohortA) { return isCohortA ? mSampleIdsCohortA : mSampleIdsCohortB; }

    public int getMaxFragmentCount() { return mMaxFragmentCount; }
    public int getMaxFragmentCount(boolean isCohortA) { return isCohortA ? mMaxFragmentCountCohortA : mMaxFragmentCountCohortB; }

    public double getAvgFragmentCount() { return mSampleIds.size() > 0 ? mTotalFragmentCount / (double)mSampleIds.size() : 0; }

    public double getAvgFragmentCount(boolean isCohortA)
    {
        return isCohortA ?
                (mSampleIdsCohortA.size() > 0 ? mTotalFragmentCountCohortA / (double)mSampleIdsCohortA.size() : 0)
                : (mSampleIdsCohortB.size() > 0 ? mTotalFragmentCountCohortB / (double)mSampleIdsCohortB.size() : 0);
    }

    public int getPositionCount(int seIndex) { return mPositionCounts[seIndex]; }
    public void addPositionCount(int seIndex, int count) { mPositionCounts[seIndex] += count; }



}
