package com.hartwig.hmftools.cobalt.calculations;

import java.util.Arrays;

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
                System.out.println("GC raw median: " + gcPail + " median: " + nextMedian);
//                MeanDepths[i] = currentMedian;//getMean(previousMedian, currentMedian, nextMedian);
                MeanDepths[i] = Arrays.stream(new double[]{previousMedian, currentMedian, nextMedian}).average().getAsDouble();
                System.out.println("Smoothed median: " + MeanDepths[i]);
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
        return gcPail != null && gcPail.mGC >= MinAllowedGc && gcPail.mGC <= MaxAllowedGc;
    }

    public double medianReadDepth(GCPail gcBucket)
    {
        if (!isAllowed(gcBucket))
        {
            return -1;
        }
        return MeanDepths[gcBucket.mGC];
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
