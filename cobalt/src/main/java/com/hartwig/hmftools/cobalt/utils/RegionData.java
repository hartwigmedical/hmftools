package com.hartwig.hmftools.cobalt.utils;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class RegionData
{
    public final int Position;
    private int mGcBucket;
    private double mMappability;

    private final List<SampleRegionData> mSampleRegionData;

    public RegionData(int position)
    {
        Position = position;
        mGcBucket = 0;
        mMappability = 0;
        mSampleRegionData = Lists.newArrayList();
    }

    public void setGcProfile(int gcBucket, double mappability)
    {
        mGcBucket = gcBucket;
        mMappability = mappability;
    }

    public void addSampleRegionData(final SampleRegionData sampleRegionData)
    {
        mSampleRegionData.add(sampleRegionData);
    }

    public String toString() { return format("%d: gcBucket(%d)", Position, mGcBucket); }
}
