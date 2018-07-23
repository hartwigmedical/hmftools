package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.min;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.calcs.DataUtils;

public class SigGroup {

    private final int mId;

    double[] mBucketRatios;
    double[] mCombinedBucketCounts;
    double mTotalCount;
    List<Integer> mSampleIds;
    double mInitialCss;
    double mWorstCss;
    double mScore;

    public SigGroup(int id, int bucketCount, double initialCSS)
    {
        mId = id;
        mBucketRatios = new double[bucketCount];
        mCombinedBucketCounts = new double[bucketCount];
        mSampleIds = Lists.newArrayList();
        mInitialCss = initialCSS;
        mWorstCss = initialCSS;
        mTotalCount = 0;
        mScore = 0;
    }

    public int id() { return mId; }
    public final double[] getBucketRatios() { return mBucketRatios; }
    public final double[] getBucketCounts() { return mCombinedBucketCounts; }
    public final List<Integer> getSampleIds() { return mSampleIds; }
    public int getTotalCount() { return (int)mTotalCount; }
    public double getInitialCss() { return mInitialCss; }
    public double getWorstCss() { return mWorstCss; }

    public double getScore() { return mScore; }
    public void setScore(double score) { mScore = score; }

    public void addSample(int sampleId, double[] bucketCounts, double css)
    {
        if(hasSampleId(sampleId))
            return;

        mSampleIds.add(sampleId);

        mWorstCss = min(mWorstCss, css);

        if(bucketCounts.length != mCombinedBucketCounts.length)
            return;

        for(int i = 0; i < mBucketRatios.length; ++i)
        {
            mCombinedBucketCounts[i] += bucketCounts[i];
        }

        calculateBucketRatios();
    }

    public boolean hasSampleId(int sampleId)
    {
        return mSampleIds.contains(sampleId);
    }

    public boolean hasSampleId(double sampleId) { return hasSampleId((int)sampleId); }

    private void calculateBucketRatios()
    {
        // calculate ratios across the buckets expressed as a percentage of the total
        // this is how Cosmic signatures are presented (but isn't necessary for NMF)
        mTotalCount = DataUtils.sumVector(mCombinedBucketCounts);

        if(mTotalCount <= 0)
            return;

        for(int i = 0; i < mBucketRatios.length; ++i)
        {
            double bucketCount = mCombinedBucketCounts[i];
            mBucketRatios[i] = bucketCount / mTotalCount;
        }
    }

}
