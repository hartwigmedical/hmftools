package com.hartwig.hmftools.data_analyser.types;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.List;

import com.google.common.collect.Lists;

public class SampleData
{
    final public int Id;

    private String mSampleName;

    private double[] mBucketCounts;
    private double[] mElevBucketCounts;
    private double[] mAllocBucketCounts;
    private double[] mUnallocBucketCounts;
    private boolean mFullyAllocated;
    private double mAllocPercent;
    private double mVarTotal;
    private double mElevatedTotal;

    private String mCancerType;

    private final List<Integer> mElevatedBuckets;
    private final List<String> mCategoryData;
    private final List<BucketGroup> mElevBucketGroups;

    public SampleData(int id)
    {
        Id = id;
        mSampleName = "";
        mCancerType = "";
        mCategoryData = Lists.newArrayList();
        mElevBucketGroups = Lists.newArrayList();
        mElevatedBuckets = Lists.newArrayList();
        mFullyAllocated = false;
        mAllocPercent = 0;
        mVarTotal = 0;
    }

    public final String getSampleName() { return mSampleName; }
    public void setSampleName(final String name) { mSampleName = name; }

    public final void setElevatedBuckets(List<Integer> buckets) { mElevatedBuckets.addAll(buckets); }
    public final List<Integer> getElevatedBuckets() { return mElevatedBuckets; }

    public final List<BucketGroup> getElevBucketGroups() { return mElevBucketGroups; }
    public void addElevBucketGroup(final BucketGroup group) { mElevBucketGroups.add(group); }

    // public final double[] getBucketCounts() { return mBucketCounts; }
    public final double[] getElevatedBucketCounts() { return mElevBucketCounts; }
    // public double getTotalCount() { return mVarTotal; }
    public double getElevatedCount() { return mElevatedTotal; }

    public final double[] getAllocBucketCounts() { return mAllocBucketCounts; }
    public final double[] getUnallocBucketCounts() { return mUnallocBucketCounts; }
    public boolean isFullyAllocated() { return mFullyAllocated; }
    public double getAllocPercent() { return mAllocPercent; }

    public void setBucketCounts(final double[] counts)
    {
        mBucketCounts = new double[counts.length];
        copyVector(counts, mBucketCounts);
        mVarTotal = sumVector(counts);
    }

    public void setElevatedBucketCounts(final double[] counts)
    {
        mElevBucketCounts = new double[counts.length];
        mAllocBucketCounts = new double[counts.length];
        mUnallocBucketCounts = new double[counts.length];
        copyVector(counts, mUnallocBucketCounts);
        copyVector(counts, mElevBucketCounts);
        mElevatedTotal = sumVector(counts);
    }

    public void allocateBucketCounts(final double[] counts)
    {
        for(int i = 0; i < counts.length; ++i)
        {
            if(mUnallocBucketCounts[i] == 0)
                continue;

            if (mUnallocBucketCounts[i] < counts[i])
            {
                mAllocBucketCounts[i] += mUnallocBucketCounts[i];
                mUnallocBucketCounts[i] = 0;
            }
            else
            {
                mUnallocBucketCounts[i] -= counts[i];
                mAllocBucketCounts[i] += counts[i];
            }
        }

        mAllocPercent = sumVector(mAllocBucketCounts)/mElevatedTotal;
    }

    public final String getCancerType() { return mCancerType; }
    public void setCancerType(final String type) { mCancerType = type; }

    public final List<String> getCategoryData() { return mCategoryData; }
    public void setCategoryData(final List<String> data) { mCategoryData.addAll(data); }


}
