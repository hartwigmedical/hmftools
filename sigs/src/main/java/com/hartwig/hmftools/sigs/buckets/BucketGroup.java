package com.hartwig.hmftools.sigs.buckets;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.sigs.DataUtils.capValue;
import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.common.sigs.DataUtils.greaterThan;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.initVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BucketGroup  {

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
    private double mMaxSimilarScore;
    private BucketGroup mMaxSimilarGroup;

    private double mPurity; // for now a percentage of sample buckets that are elevated

    public static final String BG_TYPE_BACKGROUND = "BGRD";
    public static final String BG_TYPE_MAJOR = "MAJOR";
    public static final String BG_TYPE_MINOR = "MINOR";
    public static final String BG_TYPE_UNIQUE = "UNIQUE";

    private static final Logger LOGGER = LogManager.getLogger(BucketGroup.class);

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
        mPurity = 0;

        mTotalCount = 0;
        mPotentialAllocation = 0;
        mPotentialAdjAllocation = 0;

        mCancerType = "";
        mEffects = "";
        mGroupLinks = "";
        mRefSigs = "";
        mMaxSimilarScore = 0;
        mMaxSimilarGroup = null;
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

    public boolean isBackground() { return mType.equals(BG_TYPE_BACKGROUND); }

    public int getSize() { return mBucketIds.size() * mSampleIds.size(); }

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

    public final List<Integer> getSampleIds() { return mSampleIds; }
    public final int getSampleCount() { return mSampleIds.size(); }
    public final List<Integer> getInitialSampleIds() { return mInitialSampleIds; }
    public final List<Integer> getBucketIds() { return mBucketIds; }
    public final int getBucketCount() { return mBucketIds.size(); }
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

    public void addSample(int sampleId, final double[] bucketCounts)
    {
        if(mSampleIds.contains(sampleId))
        {
            LOGGER.warn("BG({}) attempting to add sample({}) again", mId, sampleId);
            return;
        }

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

        mSampleIds.add(sampleId);

        // take a copy in case the caller retains and changes its counts
        double[] newBucketCounts = new double[bucketCounts.length];
        copyVector(bucketCounts, newBucketCounts);
        mSampleCounts.add(newBucketCounts);

        double sampleTotal = sumVector(bucketCounts);
        mSampleCountTotals.add(sampleTotal);
        mSampleCountsMap.put(sampleId, sampleTotal);
    }

    public void addSampleCounts(int samIndex, final double[] bucketCounts)
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

        setIsValid();
    }

    public final void recalcBucketRatios(double sampleWeightFactor)
    {
        // convert cumulative sample counts into ratios

        if(greaterThan(mTotalCount, 0) && mSampleIds.size() > 0)
        {
            // infer ratios from the contributing sample counts
            mBucketRatios = calcRatiosFromSampleCounts(mSampleCountTotals, mSampleCounts, sampleWeightFactor);
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
        }
    }

    public static double[] calcRatiosFromSampleCounts(final List<Double> sampleTotals, final List<double[]> sampleCountsList,
            double sampleWeightFactor)
    {
        if(sampleCountsList.isEmpty() || sampleCountsList.size() != sampleTotals.size())
        {
            return null;
        }

        int bucketCount = sampleCountsList.get(0).length;
        double[] bucketRatios = new double[bucketCount];

        double totalCount = 0;
        double[] bucketTotals = new double[bucketCount];

        for(int i = 0; i < sampleTotals.size(); ++i)
        {
            double sampleTotal = sampleTotals.get(i);
            double[] sampleCounts = sampleCountsList.get(i);

            // dampen the effect of high mutational load samples
            double sampleWF = sampleWeightFactor != 1 ? pow(sampleTotal, 1/sampleWeightFactor - 1) : 1;

            for(int j = 0; j < bucketCount; ++j)
            {
                double sbCount = sampleCounts[j] * sampleWF;
                bucketTotals[j] += sbCount;
                totalCount += sbCount;
            }
        }

        for (int i = 0; i < bucketCount; ++i)
        {
            bucketRatios[i] = bucketTotals[i] / totalCount;
        }

        return bucketRatios;
    }

    private void setIsValid()
    {
        double ratioTotal = 0;
        for (int i = 0; i < mBucketRatios.length; ++i)
        {
            if(mBucketRatios[i] < 0)
            {
                mIsValid = false;
                return;
            }

            ratioTotal += mBucketRatios[i];
        }

        if(Doubles.equal(ratioTotal, 1))
        {
            mIsValid = true;
        }
        else
        {
            mIsValid = false;
        }
    }

    public final double[] getBucketRatios()
    {
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

    public static double ratioRange(final double[] ratios, final double[] ranges, int bucket, boolean takeMin)
    {
        if(ranges == null || ranges[bucket] == 0)
            return ratios[bucket];

        double range = takeMin ? -ranges[bucket] : ranges[bucket];

        return capValue(ratios[bucket] + range, 0, 1);
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

    public double calcSampleFitScore(final SampleData sample, boolean useGross)
    {
        int samIndex = getSampleIndex(sample.Id);

        if(samIndex < 0)
            return 0;

        final double[] sampleCounts = mSampleCounts.get(samIndex);
        double sampleAlloc = mSampleCountTotals.get(samIndex);

        return calcSampleFitScore(sampleCounts, sampleAlloc, useGross);
    }

    public double calcSampleFitScore(final double[] allocCounts, double allocTotal, boolean useGross)
    {
        // the score is the ratio-weighted difference around the true ratio, with the range being the boundaries (1, -1
        // the final score is averaged by the bucket weight and so remains a score out of 1, with 0 being a perfect alignment with the ratios
        // if calculating a net score it can be negative, so between -1 and 1
        double fitScore = 0;

        if(!doublesEqual(sumVector(allocCounts), allocTotal))
            return 0;

        double weightTotal = 0;

        for(Integer bucket : mBucketIds)
        {
            double bucketWeight = allocCounts[bucket] / allocTotal;
            weightTotal += bucketWeight;

            double impliedRatio = bucketWeight;

            double ratioDiff = impliedRatio - mBucketRatios[bucket];
            double ratioRange = mBucketRatioRanges[bucket];

            if(ratioRange == 0)
                continue;

            double fitPerc = ratioDiff / ratioRange;
            fitPerc = capValue(fitPerc, -1, 1);

            if(useGross)
                fitScore += abs(fitPerc) * bucketWeight;
            else
            fitScore += fitPerc * bucketWeight;
        }

        fitScore = capValue(fitScore, -1, 1);
        fitScore /= weightTotal;

        return fitScore;
    }

    public void setMaxSimilarScore(double value) { mMaxSimilarScore = value; }
    public double getMaxSimilarScore() { return mMaxSimilarScore; }
    public final BucketGroup getMaxSimilarGroup() { return mMaxSimilarGroup; }
    public void setMaxSimilarGroup(final BucketGroup group ) { mMaxSimilarGroup = group; }

}
