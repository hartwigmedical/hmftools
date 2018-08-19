package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class BucketGroup implements Comparable<BucketGroup> {

    // keyed by a bucket pairing
    private int mId;
    private String mTag; // free-form info about the group

    private List<Integer> mSampleIds;
    private List<Integer> mBucketIds;
    private List<Integer> mInitialBucketIds;
    private List<Integer> mExtraBucketIds;

    // the bucket counts from the samples as per the specific buckets in this groiup
    private double[] mCombinedBucketCounts;
    private List<Double> mSampleCountTotals;
    private Map<Integer,Double> mSampleCountsMap;

    private boolean mBucketRatiosClean;
    private double[] mBucketRatios;
    private double[] mBucketRatioRanges;
    private List<Double> mRatioRanges;
    private double mTotalCount;

    private double mPotentialAllocation;
    private double mPotentialAdjAllocation;
    private boolean mIsSelected;

    // annotations
    private String mCancerType;
    private String mEffects;

    private double mPurity; // for now a percentage of sample buckets that are elevated
    private double mLoadFactor;
    private double mScoreOverride;

    private BucketGroup mClosestBG;
    private double mClosestBGCss;

    public BucketGroup(int id)
    {
        mId = id;
        mTag = "";

        mSampleIds = Lists.newArrayList();
        mBucketIds = Lists.newArrayList();
        mInitialBucketIds = Lists.newArrayList();
        mExtraBucketIds = Lists.newArrayList();
        mSampleCountTotals = Lists.newArrayList();
        mRatioRanges = Lists.newArrayList();
        mSampleCountsMap = new HashMap();
        mCombinedBucketCounts = null;
        mBucketRatios = null;
        mBucketRatioRanges = null;
        mBucketRatiosClean = false;
        mPurity = 0;
        mLoadFactor = 0;
        mScoreOverride = 0;

        mTotalCount = 0;
        mPotentialAllocation = 0;
        mPotentialAdjAllocation = 0;
        mIsSelected = false;

        mCancerType = "";
        mEffects = "";

        mClosestBG = null;
        mClosestBGCss = 0;
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

    public int getSize() { return mBucketIds.size() * mSampleIds.size(); }

    public double calcScore()
    {
        if(mScoreOverride > 0)
            return mScoreOverride;

        double score = sqrt(mBucketIds.size()) * mSampleIds.size();

        if(mPurity > 0)
            score *= mPurity;

        if(mLoadFactor > 0)
            score *= mLoadFactor;

        return score;
    }

    public void setScoreOverride(double score) { mScoreOverride = score; }

    public void setSelected(boolean toggle) { mIsSelected = toggle; }
    public boolean isSelected() { return mIsSelected; }

    public double getPurity() { return mPurity; }
    public void setPurity(double purity) { mPurity = purity; }
    public double getLoadFactor() { return mLoadFactor; }
    public void setLoadFactor(double value) { mLoadFactor = value; }
    public double getTotalCount() { return mTotalCount; }

    public double getAvgCount()
    {
        // per contributing sample bucket count item
        return mTotalCount/getSize();
    }

    public List<Double> getSampleCountTotals() { return mSampleCountTotals; }
    public Map<Integer, Double> getSampleCountsMap() { return mSampleCountsMap; }
    public void setSampleCountTotals(List<Double> totals) { mSampleCountTotals = totals; }

    public void setCancerType(final String type) { mCancerType = type; }
    public final String getCancerType() { return mCancerType; }

    public void setEffects(final String effects) { mEffects = effects; }
    public final String getEffects() { return mEffects; }

    public void setClosestBG(final BucketGroup bg) { mClosestBG = bg; }
    public final BucketGroup getClosestBG() { return mClosestBG; }
    public void setClosestBGCss(final double css) { mClosestBGCss = css; }
    public final double getClosestBGCss() { return mClosestBGCss; }

    public int compareTo(final BucketGroup other)
    {
        // for descending order
        return (int)round(other.calcScore() - calcScore());
    }

    public final List<Integer> getSampleIds() { return mSampleIds; }
    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final List<Integer> getInitialBucketIds() { return mInitialBucketIds; }
    public final List<Integer> getExtraBucketIds() { return mExtraBucketIds; }

    public boolean hasSample(int sampleIndex)
    {
        return mSampleIds.contains(sampleIndex);
    }

    public void clearSamples()
    {
        mSampleIds.clear();
        mSampleCountTotals.clear();
        mSampleCountsMap.clear();
        calcBucketRatios();
        mTotalCount = 0;

        for(int i = 0; i < mCombinedBucketCounts.length; ++i)
        {
            mCombinedBucketCounts[i] = 0;
        }
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

        double sampleTotal = sumVector(bucketCounts);
        mSampleCountTotals.add(sampleTotal);
        mSampleCountsMap.put(sampleId, sampleTotal);
    }

    public void merge(List<Integer> sampleIds, double[] bucketCounts)
    {
        for(Integer sample : sampleIds)
        {
            if(!mSampleIds.contains(sample))
                mSampleIds.add(sample);
        }

        for(int i = 0; i < bucketCounts.length; ++i)
        {
            // only merge counts applicable for this group's bucket collection
            if(bucketCounts[i] > 0 && mBucketIds.contains(i))
                mCombinedBucketCounts[i] += bucketCounts[i];
        }

        mBucketRatiosClean = false;
    }

    public void reduceToBucketSet(List<Integer> targetBucketIds)
    {
        int bucketIndex = 0;
        while(bucketIndex < mBucketIds.size())
        {
            int bucket = mBucketIds.get(bucketIndex);

            if(targetBucketIds.contains(bucket))
            {
                ++bucketIndex;
                continue;
            }

            // clear these out
            mCombinedBucketCounts[bucket] = 0;
            mBucketIds.remove(bucketIndex);
        }

        mBucketRatiosClean = false;
    }

    public boolean hasBucket(int bucketIndex)
    {
        return mBucketIds.contains(bucketIndex);
    }

    public boolean hasAllBuckets(final List<Integer> bucketIds)
    {
        for(Integer bucket : bucketIds)
        {
            if(!mBucketIds.contains(bucket))
                return false;
        }

        return true;
    }

    public boolean hasAnyBucket(final List<Integer> bucketIds)
    {
        for(Integer bucket : bucketIds)
        {
            if(mBucketIds.contains(bucket))
                return true;
        }

        return false;
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

    public void setBucketCounts(final double[] other)
    {
        initialise(other);
        copyVector(other, mCombinedBucketCounts);
        mBucketRatiosClean = false;
    }

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
