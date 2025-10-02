package com.hartwig.hmftools.cobalt.calculations;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;

public class GCPailsList
{
    private final List<GCPail> mBuckets = new ArrayList<>(101);

    public GCPailsList()
    {
        for (int i = 0; i < 101; i++)
        {
            mBuckets.add(new GCPail(i));
        }
    }

    public GCPail getGCPail(double gcContent)
    {
        return mBuckets.get(GCPail.bucketIndex(gcContent));
    }

    public List<GCPail> getBuckets()
    {
        return mBuckets;
    }

    public Map<GCBucket, Double> bucketToMedianReadDepth()
    {
        // todo test
        Map<GCBucket, Double> bucketToMedian = new HashMap<>();
        mBuckets.forEach(b -> bucketToMedian.put(new ImmutableGCBucket(b.mGC), b.median()));
        return bucketToMedian;
    }
}
