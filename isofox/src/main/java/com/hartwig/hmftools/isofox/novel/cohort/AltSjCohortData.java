package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

public class AltSjCohortData
{
    public final AltSpliceJunctionFile AltSJ;

    // cohort data
    private final Set<String> mSampleIds;

    private List<String> mSampleIdsCohortA;
    private List<String> mSampleIdsCohortB;
    private Map<String,List<String>> mCancerSampleIds;

    private int mTotalFragmentCount;
    private int mTotalFragmentCountCohortA;
    private int mTotalFragmentCountCohortB;

    private int mMaxFragmentCount;
    private int mMaxFragmentCountCohortA;
    private int mMaxFragmentCountCohortB;
    private final int[] mPositionCounts; // counts at the start and end

    public AltSjCohortData(final AltSpliceJunctionFile altSJ)
    {
        AltSJ = altSJ;
        mSampleIds = Sets.newHashSet();
        mSampleIdsCohortA = null;
        mSampleIdsCohortB = null;
        mCancerSampleIds = null;

        mTotalFragmentCount = 0;
        mTotalFragmentCountCohortA = 0;
        mTotalFragmentCountCohortB = 0;

        mMaxFragmentCount = 0;
        mMaxFragmentCountCohortA = 0;
        mMaxFragmentCountCohortB = 0;
        mPositionCounts = new int[SE_PAIR];
    }

    public void addSampleCount(final String sampleId, int fragCount, final String cancerType)
    {
        if(mSampleIds.contains(sampleId))
            return;

        mSampleIds.add(sampleId);
        mTotalFragmentCount += fragCount;
        mMaxFragmentCount = max(mMaxFragmentCount, fragCount);

        if(cancerType != null)
        {
            if(mCancerSampleIds == null)
                mCancerSampleIds = Maps.newHashMap();

            List<String> samples = mCancerSampleIds.get(cancerType);

            if(samples == null)
                mCancerSampleIds.put(cancerType, Lists.newArrayList(sampleId));
            else
                samples.add(sampleId);
        }
    }

    public void addSampleCount(final String sampleId, int fragCount, boolean isCohortA)
    {
        if(mSampleIdsCohortA == null)
            mSampleIdsCohortA = Lists.newArrayList();

        if(mSampleIdsCohortB == null)
            mSampleIdsCohortB = Lists.newArrayList();

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

    public final Set<String> getSampleIds() { return mSampleIds; }
    public final List<String> getCohortSampleIds(boolean isCohortA) { return isCohortA ? mSampleIdsCohortA : mSampleIdsCohortB; }
    public final Map<String,List<String>> cancerSampleIds() { return mCancerSampleIds; }

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

    public String sampleIdsStr()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        mSampleIds.forEach(x -> sj.add(x));
        return sj.toString();
    }


}
