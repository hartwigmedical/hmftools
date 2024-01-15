package com.hartwig.hmftools.common.genome.gc;

import java.util.Map;

public class GCMedianReadDepth
{
    private final double mMean;
    private final double mMedian;
    private final Map<GCBucket, Double> mMedianReadDepthPerGCBucket;

    public GCMedianReadDepth(final double mean, final double median, final Map<GCBucket,Double> medianReadCountPerGCBucket)
    {
        mMean = mean;
        mMedian = median;
        mMedianReadDepthPerGCBucket = medianReadCountPerGCBucket;
    }

    public double meanReadDepth()
    {
        return mMean;
    }

    public double medianReadDepth()
    {
        return mMedian;
    }

    public double medianReadDepth(final GCBucket bucket)
    {
        return mMedianReadDepthPerGCBucket.getOrDefault(bucket, -1.0);
    }

    public double medianReadDepth(final GCProfile profile)
    {
        return medianReadDepth(GCBucket.create(profile));
    }
}
