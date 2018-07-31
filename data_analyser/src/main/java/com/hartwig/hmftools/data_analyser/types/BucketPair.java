package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import java.util.List;

import com.google.common.collect.Lists;

public class BucketPair {

    // keyed by a bucket pairing
    int mBpId;
    int mBucketA;
    int mBucketB;

    // a list of all samples whose ratio is within a similar range
    List<Integer> mSampleIds;
    List<Double> mBPRatios;

    double mHighRatio;
    double mLowRatio;
    double mAverageRatio;

    public static double RATIO_MAX = 100;
    public static double RATIO_MIN = 1 / RATIO_MAX;
    private static double RATIO_MARGIN = 0.05; // 5%, used in comparisons and rounding

    public BucketPair(int bpId, int bucketA, int bucketB)
    {
        mBpId = bpId;
        mBucketA = bucketA;
        mBucketB = bucketB;

        mSampleIds = Lists.newArrayList();
        mBPRatios = Lists.newArrayList();

        mLowRatio = RATIO_MAX;
        mHighRatio = RATIO_MIN;
        mAverageRatio = 0;
    }

    public int getId() { return mBpId; }
    public int getBucketA() { return mBucketA; }
    public int getBucketB() { return mBucketB; }

    public final List<Integer> getSampleIds() { return mSampleIds; }
    public final List<Double> getRatios() { return mBPRatios; }

    public void addSampleRatio(int sampleIndex, double bpRatio)
    {
        mSampleIds.add(sampleIndex);
        mBPRatios.add(bpRatio);

        mLowRatio = Math.min(mLowRatio, bpRatio);
        mHighRatio = Math.max(mHighRatio, bpRatio);
        mAverageRatio = (mAverageRatio * (mBPRatios.size()-1) + bpRatio) / mBPRatios.size();
    }

    public double getHighRatio() { return mHighRatio; }
    public double getLowRatio() { return mLowRatio; }

    public double getAverageRatio() { return mAverageRatio; }

    public double calcRangePercent()
    {
        double range = mHighRatio - mLowRatio;

        if(range <= 0)
            return 0;

        return range / mAverageRatio;
    }

    public double calcAverageRatio()
    {
        if(mBPRatios.size() == 0)
            return 0;

        double total = 0;
        for(final Double ratio : mBPRatios)
        {
            total += ratio;
        }

        return total / mBPRatios.size();
    }

    public boolean isWithinRange(double bpRatio)
    {
        if(bpRatio < mAverageRatio * (1-RATIO_MARGIN))
            return false;
        else if(bpRatio > mAverageRatio * (1+RATIO_MARGIN))
            return false;
        else
            return true;
    }

    public boolean hasSameBucket(final BucketPair other)
    {
        return other.getBucketA() == mBucketA
            || other.getBucketA() == mBucketB
            || other.getBucketB() == mBucketA
            || other.getBucketB() == mBucketB;
    }

    public List<Integer> getSharedSamples(final List<Integer> otherSamples)
    {
        List<Integer> samplesList = Lists.newArrayList();

        // early exit if the boundaries don't overlap, since the sampleIds are in ascending order
        if(mSampleIds.get(0) > otherSamples.get(otherSamples.size()-1)
        || mSampleIds.get(mSampleIds.size()-1) < otherSamples.get(0))
        {
            return samplesList;
        }

        for(final Integer sample : otherSamples)
        {
            if(hasSample(sample))
                samplesList.add(sample);
        }

        return samplesList;
    }

    public boolean hasSample(int sampleIndex)
    {
        return mSampleIds.contains(sampleIndex);
    }

    private static double roundRatio(double ratio)
    {
        // scaled round mechanism
        double log10 = floor(log10(ratio * RATIO_MARGIN));
        double roundFactor = pow(10,log10);
        double ratioRounded = round(ratio/roundFactor)*roundFactor;
        return ratioRounded;
    }

}
