package com.hartwig.hmftools.data_analyser.types;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;

import java.util.List;

import com.google.common.collect.Lists;

public class SampleData
{
    final public int Id;

    private String mSampleName;

    private double[] mBucketCounts;

    private String mCancerType;

    final List<String> mCategoryData;

    public SampleData(int id)
    {
        Id = id;
        mSampleName = "";
        mCancerType = "";
        mCategoryData = Lists.newArrayList();
    }

    public final String getSampleName() { return mSampleName; }
    public void setSampleName(final String name) { mSampleName = name; }

    public final double[] getBucketCounts() { return mBucketCounts; }
    public void setBucketCounts(final double[] counts)
    {
        mBucketCounts = new double[counts.length];
        copyVector(counts, mBucketCounts);
    }

    public final String getCancerType() { return mCancerType; }
    public void setCancerType(final String type) { mCancerType = type; }

    public final List<String> getCategoryData() { return mCategoryData; }
    public void setCategoryData(final List<String> data) { mCategoryData.addAll(data); }


}
