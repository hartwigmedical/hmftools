package com.hartwig.hmftools.common.genome.gc;

import java.util.Map;

public class GCMedianReadCount
{
    private final double mMean;
    private final double mMedian;
    private final Map<GCBucket, Double> mMedianReadCountPerGCBucket;

    public GCMedianReadCount(final double mean, final double median, final Map<GCBucket,Double> medianReadCountPerGCBucket)
    {
        mMean = mean;
        mMedian = median;
        mMedianReadCountPerGCBucket = medianReadCountPerGCBucket;
    }

    public double meanReadCount()
    {
        return mMean;
    }

    public double medianReadCount()
    {
        return mMedian;
    }

    public double medianReadCount(final GCBucket bucket)
    {
        return mMedianReadCountPerGCBucket.getOrDefault(bucket, -1.0);
    }

    public double medianReadCount(final GCProfile profile)
    {
        return medianReadCount(GCBucket.create(profile));
    }
}
