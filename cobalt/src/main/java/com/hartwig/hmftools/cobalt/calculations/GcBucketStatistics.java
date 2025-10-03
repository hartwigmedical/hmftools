package com.hartwig.hmftools.cobalt.calculations;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;

public class GcBucketStatistics
{
    private final int MinAllowedGc;
    private final int MaxAllowedGc;
    private final double[] MeanDepths = new double[101];

    public GcBucketStatistics(GCPailsList gcPailsList, int minAllowedGc, int maxAllowedGc)
    {
        Preconditions.checkArgument(minAllowedGc > 0);
        Preconditions.checkArgument(maxAllowedGc < 100);
        this.MinAllowedGc = minAllowedGc;
        this.MaxAllowedGc = maxAllowedGc;
        double[] window = {-1.0, -1.0, -1.0};
        for(int i = 0; i < MeanDepths.length - 1; i++)
        {
            GCPail gcPail = gcPailsList.getBuckets().get(i+1);
            window[0] = window[1];
            window[1] = window[2];
            window[2] = gcPail.median();
            if(isAllowed(gcPail))
            {
                MeanDepths[i] = Arrays.stream(window).average().orElse(-1.0);
            }
            else
            {
                MeanDepths[i] = -1.0;
            }
        }
        MeanDepths[100] = -1.0;
    }

    public Map<GCBucket, Double> bucketToMedianReadDepth()
    {
        // todo test
        Map<GCBucket, Double> bucketToMedian = new HashMap<>();
        for (int i=0; i<101; i++)
        {
            if(MeanDepths[i] > 0)
            {
                bucketToMedian.put(new ImmutableGCBucket(i), MeanDepths[i]);
            }
        }
        return bucketToMedian;
    }

    public boolean isAllowed(GCPail gcPail)
    {
        return gcPail != null && gcPail.mGC > MinAllowedGc && gcPail.mGC <= MaxAllowedGc;
    }

    public double medianReadDepth(GCPail gcBucket)
    {
        if (!isAllowed(gcBucket))
        {
            return -1;
        }
        return MeanDepths[gcBucket.mGC];
    }
}
