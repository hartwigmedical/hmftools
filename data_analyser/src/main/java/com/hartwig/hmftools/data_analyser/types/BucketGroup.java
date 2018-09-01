package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVectors;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class BucketGroup implements Comparable<BucketGroup> {

    // keyed by a bucket pairing
    private int mId;
    private String mTag; // free-form info about the group

    private List<Integer> mSampleIds;
    private List<Integer> mInitialSampleIds; // those which led to the creation of the group
    private List<Integer> mBucketIds;
    private List<Integer> mInitialBucketIds;
    private List<Integer> mExtraBucketIds;

    // the bucket counts from the samples as per the specific buckets in this groiup
    private double[] mCombinedBucketCounts;
    private List<Double> mSampleCountTotals;
    private Map<Integer,Double> mSampleCountsMap;
    private List<double[]> mSampleCounts;

    private boolean mBucketRatiosClean;
    private double[] mBucketRatios;
    private double[] mBucketRatioRanges;
    private List<Double> mRatioRanges;
    private double mTotalCount;

    private double mPotentialAllocation;
    private double mPotentialAdjAllocation;

    // annotations
    private String mCancerType;
    private String mEffects;

    private double mPurity; // for now a percentage of sample buckets that are elevated
    private double mScoreOverride;

    public static String BG_BACKGROUND_TYPE = "background";

    public BucketGroup(int id)
    {
        mId = id;
        mTag = "";

        mSampleIds = Lists.newArrayList();
        mInitialSampleIds = Lists.newArrayList();
        mBucketIds = Lists.newArrayList();
        mInitialBucketIds = Lists.newArrayList();
        mExtraBucketIds = Lists.newArrayList();
        mSampleCountTotals = Lists.newArrayList();
        mSampleCounts = Lists.newArrayList();
        mRatioRanges = Lists.newArrayList();
        mSampleCountsMap = new HashMap();
        mCombinedBucketCounts = null;
        mBucketRatios = null;
        mBucketRatioRanges = null;
        mBucketRatiosClean = false;
        mPurity = 0;
        mScoreOverride = 0;

        mTotalCount = 0;
        mPotentialAllocation = 0;
        mPotentialAdjAllocation = 0;

        mCancerType = "";
        mEffects = "";
    }

    private void initialise(final double[] counts)
    {
        if(mBucketRatios == null)
            mBucketRatios = new double[counts.length];

        if(mBucketRatioRanges == null)
            mBucketRatioRanges = new double[counts.length];

        if(mCombinedBucketCounts == null)
            mCombinedBucketCounts = new double[counts.length];
    }

    public int getId() { return mId; }

    public String getTag() { return mTag; }
    public void setTag(final String tag) { mTag = tag; }
    public boolean isBackground() { return mTag.equals(BG_BACKGROUND_TYPE); }

    public int getSize() { return mBucketIds.size() * mSampleIds.size(); }

    public double calcScore()
    {
        if(mScoreOverride > 0)
            return mScoreOverride;

        double score = sqrt(mBucketIds.size()) * mSampleIds.size();

        if(mPurity > 0)
            score *= mPurity;

        return score;
    }

    public double getPurity() { return mPurity; }
    public void setPurity(double purity) { mPurity = purity; }
    public double getTotalCount() { return mTotalCount; }

    public double getAvgCount()
    {
        // per contributing sample bucket count item
        return mTotalCount/getSize();
    }

    public List<Double> getSampleCountTotals() { return mSampleCountTotals; }
    public Map<Integer, Double> getSampleCountsMap() { return mSampleCountsMap; }
    public List<double[]> getSampleCounts() { return mSampleCounts; }
    public double getSampleCount(Integer sampleId)
    {
        if(mSampleCountsMap.containsKey(sampleId))
            return mSampleCountsMap.get(sampleId);
        else
            return 0;
    }

    public void setCancerType(final String type) { mCancerType = type; }
    public final String getCancerType() { return mCancerType; }

    public void setEffects(final String effects) { mEffects = effects; }
    public final String getEffects() { return mEffects; }

    public int compareTo(final BucketGroup other)
    {
        // for descending order
        return (int)round(other.calcScore() - calcScore());
    }

    public final List<Integer> getSampleIds() { return mSampleIds; }
    public final List<Integer> getInitialSampleIds() { return mInitialSampleIds; }
    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final List<Integer> getInitialBucketIds() { return mInitialBucketIds; }
    public final List<Integer> getExtraBucketIds() { return mExtraBucketIds; }

    public boolean hasSample(Integer sampleId)
    {
        return mSampleIds.contains(sampleId);
    }

    public boolean hasInitialSample(Integer sampleId)
    {
        return mInitialSampleIds.contains(sampleId);
    }

    public void clearSamples()
    {
        mSampleIds.clear();
        mSampleCountTotals.clear();
        mSampleCountsMap.clear();
        mSampleCounts.clear();
        calcBucketRatios();
        mTotalCount = 0;

        initVector(mCombinedBucketCounts, 0);
    }

    public void addInitialSample(int sampleId)
    {
        if(mInitialSampleIds.contains(sampleId))
            return;

        mInitialSampleIds.add(sampleId);
    }

    public void addSample(int sampleId, double[] bucketCounts)
    {
        addSample(sampleId, bucketCounts, true);
    }

    public void addSample(int sampleId, double[] bucketCounts, boolean reqRatioRecalc)
    {
        if(mSampleIds.contains(sampleId))
            return;

        initialise(bucketCounts);

        for(Integer bucketId : mBucketIds)
        {
            mCombinedBucketCounts[bucketId] += bucketCounts[bucketId];
            mTotalCount += bucketCounts[bucketId];
        }

        if(reqRatioRecalc)
            mBucketRatiosClean = false;

        mSampleIds.add(sampleId);
        mSampleCounts.add(bucketCounts);

        double sampleTotal = sumVector(bucketCounts);
        mSampleCountTotals.add(sampleTotal);
        mSampleCountsMap.put(sampleId, sampleTotal);
    }

    public boolean removeSampleAllocation(final SampleData sample, boolean removePotentialAlloc)
    {
        int samIndex = -1;
        for(int index = 0; index < mSampleIds.size(); ++index)
        {
            if(mSampleIds.get(index) == sample.Id)
            {
                samIndex = index;
                break;
            }
        }

        if(samIndex == -1)
            return false;

        if(removePotentialAlloc)
        {
            double sampleAlloc = mSampleCountTotals.get(samIndex);
            mPotentialAllocation -= sampleAlloc;
            mPotentialAdjAllocation -= sampleAlloc / sample.getElevatedCount();
        }

        final double[] sampleCounts = mSampleCounts.get(samIndex);
        for(Integer bucketId : mBucketIds)
        {
            mCombinedBucketCounts[bucketId] = sampleCounts[bucketId];
            mTotalCount += sampleCounts[bucketId];
        }

        mSampleIds.remove(samIndex);
        mSampleCountTotals.remove(samIndex);
        mSampleCountsMap.remove(sample.Id);
        mSampleCounts.remove(samIndex);

        if(mSampleIds.size() != mSampleCountTotals.size() || mSampleIds.size() != mSampleCountsMap.size() || mSampleIds.size() != mSampleCounts.size())
        {
            return false;
        }

        return true;
    }

    public boolean hasBucket(int bucketIndex)
    {
        return mBucketIds.contains(bucketIndex);
    }

    public void addBuckets(List<Integer> bucketIds)
    {
        for(Integer bucket : bucketIds)
        {
            addBucket(bucket,true);
        }
    }

    public void addBucket(int bucketId, boolean isInitial)
    {
        if(mBucketIds.contains(bucketId))
            return;

        mBucketIds.add(bucketId);

        if(isInitial)
            mInitialBucketIds.add(bucketId);
        else
            mExtraBucketIds.add(bucketId);
    }

    public final double[] getBucketCounts() { return mCombinedBucketCounts; }

    public void setBucketRatios(final double[] other)
    {
        initialise(other);

        copyVector(other, mBucketRatios);
        mTotalCount = sumVector(mCombinedBucketCounts);
        mBucketRatiosClean = true;
    }

    public final void calcBucketRatios()
    {
        if(mBucketRatiosClean)
            return;

        mTotalCount = sumVector(mCombinedBucketCounts);

        for (int i = 0; i < mBucketRatios.length; ++i)
        {
            mBucketRatios[i] = mCombinedBucketCounts[i] / mTotalCount;
        }

        double ratioTotal = sumVector(mBucketRatios);
        if(doublesEqual(ratioTotal, 1))
            mBucketRatiosClean = true;
    }

    public final double[] getBucketRatios()
    {
        calcBucketRatios();
        return mBucketRatios;
    }

    // public final double[] getBucketRatioRanges() { return mBucketRatioRanges; }
    public final List<Double> getBucketRatioRanges() { return mRatioRanges; }

    public void setBucketRatioRanges(final List<Double> ranges)
    {
        mRatioRanges.clear();
        mRatioRanges.addAll(ranges);

        if (mBucketIds.size() != ranges.size())
            return;

        for (int bIndex = 0; bIndex < mBucketIds.size(); ++bIndex)
        {
            Integer bucket = mBucketIds.get(bIndex);
            mBucketRatioRanges[bucket] = ranges.get(bIndex);
        }
    }

    public double getPotentialAllocation() { return mPotentialAllocation; }
    public void addPotentialAllocation(double count) { mPotentialAllocation += count; }
    public void resetPotentialAllocation()
    {
        mPotentialAllocation = 0;
        mPotentialAdjAllocation = 0;
    }

    public double getPotentialAdjAllocation() { return mPotentialAdjAllocation; }
    public void addPotentialAdjAllocation(double count) { mPotentialAdjAllocation += count; }
}
