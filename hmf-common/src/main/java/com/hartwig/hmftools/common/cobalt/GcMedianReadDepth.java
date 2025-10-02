package com.hartwig.hmftools.common.cobalt;

import java.util.Map;

import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.GCProfile;

public class GcMedianReadDepth
{
    private final double mMean;
    private final double mMedian;
    private final Map<GCBucket, Double> mMedianReadDepthPerGCBucket;

    public GcMedianReadDepth(final double mean, final double median, final Map<GCBucket,Double> medianReadCountPerGCBucket)
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
}
