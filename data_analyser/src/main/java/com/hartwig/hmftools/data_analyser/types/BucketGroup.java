package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.greaterThan;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;

public class BucketGroup implements Comparable<BucketGroup> {

    // keyed by a bucket pairing
    private int mId;
    private String mTag; // free-form info about the group
    private String mType;
    private boolean mIsValid;

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
    private double[] mInitialBucketRatios;
    private double[] mBucketRatioRanges;
    private double mTotalCount;

    private double mPotentialAllocation;
    private double mPotentialAdjAllocation;

    // annotations
    private String mCancerType;
    private String mEffects;
    private String mGroupLinks;
    private String mRefSigs;

    private double mPurity; // for now a percentage of sample buckets that are elevated
    private double mScoreOverride;

    public static String BG_TYPE_BACKGROUND = "BGRD";
    public static String BG_TYPE_MAJOR = "MAJOR";
    public static String BG_TYPE_MINOR = "MINOR";
    public static String BG_TYPE_UNIQUE = "UNIQUE";

    public BucketGroup(int id)
    {
        mId = id;
        mTag = "";
        mType = "";
        mIsValid = false;

        mSampleIds = Lists.newArrayList();
        mInitialSampleIds = Lists.newArrayList();
        mBucketIds = Lists.newArrayList();
        mInitialBucketIds = Lists.newArrayList();
        mExtraBucketIds = Lists.newArrayList();
        mSampleCountTotals = Lists.newArrayList();
        mSampleCounts = Lists.newArrayList();
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
        mGroupLinks = "";
        mRefSigs = "";
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
    public boolean isValid() { return mIsValid; }

    public String getTag() { return mTag; }
    public void setTag(final String tag) { mTag = tag; }

    public String getGroupType() { return mType; }
    public void setGroupType(final String type) { mType = type; }

    public boolean isBackground() { return mTag.equals(BG_TYPE_BACKGROUND); }

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
    public double getTotalCount() { return mTotalCount; }

    public double getAvgCount()
    {
        // per contributing sample bucket count item
        return mTotalCount/getSize();
    }

    public List<Double> getSampleCountTotals() { return mSampleCountTotals; }
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

    public void setEffects(final String value) { mEffects = value; }
    public final String getEffects() { return mEffects; }

    public void addGroupLinks(final String value)
    {
        if(!mGroupLinks.isEmpty())
            mGroupLinks += ";";

        mGroupLinks += value;
    }

    public final String getGroupLinks() { return mGroupLinks; }

    public void addRefSig(final String value)
    {
        if(!mRefSigs .isEmpty())
            mRefSigs  += ";";

        mRefSigs += value;
    }

    public final String getRefSigs() { return mRefSigs; }

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

    public void clearSamples()
    {
        mSampleIds.clear();
        mSampleCountTotals.clear();
        mSampleCountsMap.clear();
        mSampleCounts.clear();
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
            if(bucketCounts[bucketId] < 0)
            {
                mIsValid = false;
                return;
            }

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

    public void addSampleCounts(int samIndex, double[] bucketCounts)
    {
        // add to existing counts for an existing sample
        if(samIndex < 0 || samIndex >= mSampleIds.size())
            return;

        int sampleId = mSampleIds.get(samIndex);

        double[] existingCounts = mSampleCounts.get(samIndex);

        double countsTotal = 0;
        for(Integer bucketId : mBucketIds)
        {
            mCombinedBucketCounts[bucketId] += max(bucketCounts[bucketId], 0);
            existingCounts[bucketId] += bucketCounts[bucketId];
            countsTotal += bucketCounts[bucketId];
        }

        mTotalCount += countsTotal;
        double newTotal = mSampleCountTotals.get(samIndex) + countsTotal;
        mSampleCountTotals.set(samIndex, newTotal);
        mSampleCountsMap.put(sampleId, newTotal);
    }

    public int getSampleIndex(int sampleId)
    {
        for (int index = 0; index < mSampleIds.size(); ++index)
        {
            if (mSampleIds.get(index) == sampleId)
            {
                return index;
            }
        }

        return -1;
    }

    public boolean removeSampleAllocation(final SampleData sample, int samIndex, boolean removePotentialAlloc)
    {
        if(samIndex == -1)
        {
            samIndex = getSampleIndex(sample.Id);

            if (samIndex == -1)
                return false;
        }

        if(removePotentialAlloc)
        {
            double sampleAlloc = mSampleCountTotals.get(samIndex);
            mPotentialAllocation -= sampleAlloc;
            mPotentialAdjAllocation -= sampleAlloc * (sampleAlloc / sample.getElevatedCount());

            mBucketRatiosClean = false;
        }

        final double[] sampleCounts = mSampleCounts.get(samIndex);
        for(Integer bucketId : mBucketIds)
        {
            mCombinedBucketCounts[bucketId] = max(mCombinedBucketCounts[bucketId] - sampleCounts[bucketId], 0);
            mTotalCount = max(mTotalCount - sampleCounts[bucketId], 0);
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

        setIsValid();

        if(mIsValid)
            mBucketRatiosClean = true;
    }

    public final void recalcBucketRatios()
    {
        calcBucketRatios(true);
    }

    private final void calcBucketRatios(boolean forceRecalc)
    {
        if(mBucketRatiosClean && !forceRecalc)
            return;

        if(greaterThan(mTotalCount, 0) && mSampleIds.size() > 0)
        {
            // infer ratios from the contributing sample counts
            mTotalCount = sumVector(mCombinedBucketCounts);

            for (int i = 0; i < mBucketRatios.length; ++i)
            {
                mBucketRatios[i] = mCombinedBucketCounts[i] / mTotalCount;
            }
        }
        else if(mInitialBucketRatios != null)
        {
            copyVector(mInitialBucketRatios, mBucketRatios);
        }

        setIsValid();

        if(mIsValid)
        {
            if(mInitialBucketRatios == null)
            {
                mInitialBucketRatios = new double[mBucketRatios.length];
                copyVector(mBucketRatios, mInitialBucketRatios);
            }

            mBucketRatiosClean = true;
        }
    }

    private void setIsValid()
    {
        double ratioTotal = 0;
        for (int i = 0; i < mBucketRatios.length; ++i)
        {
            if(mBucketRatios[i] < 0)
            {
                mBucketRatiosClean = false;
                mIsValid = false;
                return;
            }

            ratioTotal += mBucketRatios[i];
        }

        if(Doubles.equal(ratioTotal, 1))
        {
            mBucketRatiosClean = true;
            mIsValid = true;
        }
        else
        {
            mIsValid = false;
        }
    }

    public final double[] getBucketRatios()
    {
        calcBucketRatios(false);
        return mBucketRatios;
    }

    public final double[] getRatioRanges() { return mBucketRatioRanges; }

    public void setBucketRatioRanges(final double[] ratioRanges)
    {
        if(ratioRanges != null)
            copyVector(ratioRanges, mBucketRatioRanges);
    }

    public void setRatioRangePerc(double rangePerc)
    {
        for(Integer bucket : mBucketIds)
        {
            mBucketRatioRanges[bucket] = mBucketRatios[bucket] * rangePerc;
        }
    }

    public static double ratioRange(final double[] ranges, int bucket, boolean takeMin)
    {
        if(ranges == null || ranges[bucket] == 0)
            return 0;

        return takeMin ? -ranges[bucket] : ranges[bucket];
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
