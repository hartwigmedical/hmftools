package com.hartwig.hmftools.cobalt.calculations;

import com.google.common.base.Preconditions;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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
        double previousMedian;
        double currentMedian = -1.0;
        double nextMedian = -1.0;
        for(int i = 0; i < MeanDepths.length - 1; i++)
        {
            GCPail gcPail = gcPailsList.getBuckets().get(i+1);
            previousMedian = currentMedian;
            currentMedian = nextMedian;
            nextMedian = gcPail.median();
            if(isAllowed(gcPail))
            {
                MeanDepths[i] = getMean(previousMedian, currentMedian, nextMedian);
            }
            else
            {
                MeanDepths[i] = -1.0;
            }
        }
        MeanDepths[100] = -1.0;
    }

    public boolean isAllowed(GCPail gcPail)
    {
        return gcPail.mGC >= MinAllowedGc && gcPail.mGC <= MaxAllowedGc;
    }

    public double medianReadDepth(int gcBucket)
    {
        Preconditions.checkArgument(gcBucket >= 0);
        Preconditions.checkArgument(gcBucket <= 100);
        return MeanDepths[gcBucket];
    }

    public double medianReadDepthAcrossInRangeBuckets()
    {
        return 0;
    }

    public double meanReadDepthAcrossInRangeBuckets()
    {
        return 0;
    }

    private double getMean(double previous, double current, double next)
    {
        if(current < 0)
        {
            return 0.0;
        }
        DescriptiveStatistics stats = new DescriptiveStatistics();
        if(previous >= 0)
        {
            stats.addValue(previous);
        }
        stats.addValue(current);
        if(next >= 0)
        {
            stats.addValue(next);
        }
        return stats.getMean();
    }
}
