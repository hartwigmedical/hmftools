package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.round;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.List;

import com.google.common.collect.Lists;

public class BucketGroup implements Comparable<BucketGroup> {

    // keyed by a bucket pairing
    int mId;

    List<Integer> mSampleIds;
    List<Integer> mBucketIds;
    List<Integer> mInitialBucketIds;
    List<Integer> mExtraBucketIds;

    // the bucket counts from the samples as per the specific buckets in this groiup
    double[] mCombinedBucketCounts;

    boolean mBucketRatiosClean;
    double[] mBucketRatios;

    // annotations
    private String mCancerType;
    private String mEffects;

    double mPurity; // for now a percentage of sample buckets that are elevated

    public BucketGroup(int id)
    {
        mId = id;

        mSampleIds = Lists.newArrayList();
        mBucketIds = Lists.newArrayList();
        mInitialBucketIds = Lists.newArrayList();
        mExtraBucketIds = Lists.newArrayList();
        mCombinedBucketCounts = null;
        mBucketRatios = null;
        mBucketRatiosClean = false;
        mPurity = 0;

        mCancerType = "Unclear";
        mEffects = "";
    }

    public int getId() { return mId; }

    public int getSize() { return mBucketIds.size() * mSampleIds.size(); }

    public double calcScore()
    {
        if(mPurity == 0)
            return getSize();

        return getSize() * mPurity;
    }

    public double getPurity() { return mPurity; }
    public void setPurity(double purity) { mPurity = purity; }

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
    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final List<Integer> getInitialBucketIds() { return mInitialBucketIds; }
    public final List<Integer> getExtraBucketIds() { return mExtraBucketIds; }

    public boolean hasSample(int sampleIndex)
    {
        return mSampleIds.contains(sampleIndex);
    }

    public void addSample(int sampleId, double[] bucketCounts)
    {
        if(mSampleIds.contains(sampleId))
            return;

        if(mCombinedBucketCounts == null)
            mCombinedBucketCounts = new double[bucketCounts.length];
        else if(bucketCounts.length != mCombinedBucketCounts.length)
            return;

        for(int i = 0; i < bucketCounts.length; ++i)
        {
            mCombinedBucketCounts[i] += bucketCounts[i];
        }

        mBucketRatiosClean = false;

        mSampleIds.add(sampleId);
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

    public void addBuckets(List<Integer> bucketIds)
    {
        for(Integer bucket : bucketIds)
        {
            addBucket(bucket,true);
        }
    }

    public void addBucket(int bucketId, double[] bucketCounts, boolean isInitial)
    {
        if(mBucketIds.contains(bucketId))
            return;

        addBucket(bucketId, isInitial);

        mCombinedBucketCounts[bucketId] = sumVector(bucketCounts);
    }

    private void addBucket(int bucketId, boolean isInitial)
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
        copyVector(other, mCombinedBucketCounts);
        mBucketRatiosClean = false;
    }

    public final double[] getBucketRatios()
    {
        if(mBucketRatios == null)
            mBucketRatios = new double[mCombinedBucketCounts.length];

        if(!mBucketRatiosClean)
        {
            double totalCount = sumVector(mCombinedBucketCounts);

            for (int i = 0; i < mBucketRatios.length; ++i)
            {
                mBucketRatios[i] = mCombinedBucketCounts[i] / totalCount;
            }

            double ratioTotal = sumVector(mBucketRatios);
            if(doublesEqual(ratioTotal, 1))
                mBucketRatiosClean = true;
        }

        return mBucketRatios;
    }

    public static List<Integer> getMatchingBucketList(final List<Integer> bl1, final List<Integer> bl2)
    {
        // gets union/common set
        List<Integer> matchedList = Lists.newArrayList();

        for(Integer bucket : bl1)
        {
            if(bl2.contains(bucket))
                matchedList.add(bucket);
        }

        return matchedList;
    }

    public static List<Integer> getDiffBucketList(final List<Integer> bl1, final List<Integer> bl2)
    {
        // returns list of buckets in 1 but not in 2
        List<Integer> diffList = Lists.newArrayList();

        for(Integer bucket : bl1)
        {
            if(!bl2.contains(bucket))
                diffList.add(bucket);
        }

        return diffList;
    }

    public static List<Integer> getCombinedBuckets(final List<Integer> bl1, final List<Integer> bl2)
    {
        // gets super set, including non-common buckets
        List<Integer> combinedSet = Lists.newArrayList();

        for(Integer bucket : bl1)
        {
            if(!combinedSet.contains(bucket))
                combinedSet.add(bucket);
        }

        for(Integer bucket : bl2)
        {
            if(!combinedSet.contains(bucket))
                combinedSet.add(bucket);
        }

        return combinedSet;
    }

}
