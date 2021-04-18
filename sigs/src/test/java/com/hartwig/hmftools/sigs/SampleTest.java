package com.hartwig.hmftools.sigs;

import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MAX_NOISE_ALLOC_PERCENT;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sigs.buckets.SampleData;

import org.junit.Test;

public class SampleTest
{
    @Test
    public void testPotentialAllocations()
    {
        MAX_NOISE_ALLOC_PERCENT = 0.2; // max noise of 20% of the total var count

        SampleData sample = new SampleData(0);

        int bucketCount = 4;
        double[] counts = {10, 20, 30, 40};
        double[] elevCounts = counts;

        double[] noiseCounts = new double[bucketCount];

        sample.setBucketCounts(counts);
        sample.setElevatedBucketCounts(elevCounts, noiseCounts);

        // first test perfect allocation, with no noise or ratio ranges
        double[] bucketRatios = {0.1, 0.2, 0.3, 0.4};
        double[] ratioRanges = new double[bucketCount];

        // just first 2 buckets
        List<Integer> requiredBuckets = Lists.newArrayList();
        requiredBuckets.add(0);
        requiredBuckets.add(1);

        double[] allocCounts = new double[bucketCount];
        double allocTotal = sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(30.0, sumVector(allocCounts));

        requiredBuckets.add(2);
        requiredBuckets.add(3);

        allocTotal = sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(100.0, sumVector(allocCounts));

        // now test with noise in each bucket

        double noisePerBucket = 10;
        for(int i = 0;i < bucketCount; ++i)
        {
            noiseCounts[i] = noisePerBucket;
        }

        sample.setElevatedBucketCounts(elevCounts, noiseCounts);

        sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(100.0, sumVector(allocCounts)); // perfect allocation skips boosting allocs with noise

        // test again with one ratio being the limiting factor
        for(int i = 0;i < bucketCount; ++i)
        {
            bucketRatios[i] = 0.25;
        }

        // bucket 1 limits the others (4 x 10 out of 100), so noise is also limited to 80% of the max 20 available
        allocTotal = sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(20.0 + 3 * 20, allocTotal);
        assertEquals(20.0, allocCounts[0]);

        // test again with 2 buckets limited and competing for the same noise
        bucketRatios[0] = 0.2;
        bucketRatios[1] = 0.4;
        bucketRatios[2] = 0.2;
        bucketRatios[3] = 0.2;

        allocTotal = sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(15.0 + 20 + 10 + 2 * 15, allocTotal); // all of buckets 1 & 2, full noise, plus ratio-limited for 3 & 5
        assertEquals(10.0 + 1/3.0 * 15.0, allocCounts[0], 0.001);
        assertEquals(20.0 + 2/3.0 * 15.0, allocCounts[1], 0.001);

        // finally test with ratio ranges
        for (int i = 0; i < bucketCount; ++i)
        {
            ratioRanges[i] = 0.05;
        }

        bucketRatios[0] = 0.06;
        bucketRatios[1] = 0.15;
        bucketRatios[2] = 0.35;
        bucketRatios[3] = 0.38;

        allocTotal = sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(100.0, allocTotal);
    }

    @Test
    public void testAllocations()
    {
        MAX_NOISE_ALLOC_PERCENT = 0.2; // max noise of 20% of the total var count

        SampleData sample = new SampleData(0);

        int bucketCount = 4;
        double[] counts = { 100, 200, 300, 400 };
        double[] elevCounts = counts;
        double[] noiseCounts = new double[bucketCount];

        double noisePerBucket = 20;
        List<Integer> requiredBuckets = Lists.newArrayList();

        for (int i = 0; i < bucketCount; ++i)
        {
            requiredBuckets.add(i);
            noiseCounts[i] = noisePerBucket;
        }

        sample.setBucketCounts(counts);
        sample.setElevatedBucketCounts(elevCounts, noiseCounts);

        // first test perfect allocation, with no noise or ratio ranges
        double[] bucketRatios = { 0.1, 0.2, 0.3, 0.4 };
        double[] ratioRanges = new double[bucketCount];

        double[] allocCounts = new double[bucketCount];
                sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(1000.0, sumVector(allocCounts));

        double allocTotal = sample.allocateBucketCounts(allocCounts, 0.5);
        assertEquals(1000.0, allocTotal);
        assertEquals(1000.0, sample.getAllocatedCount());
        assertEquals(1.00, sample.getAllocPercent());
        assertEquals(0.0, sample.getAllocNoise());

        // remove allocation then try again with ratio and ranges to cover counts exactly
        sample.clearAllocations(true);

        // now test with imperfect ratios and use of noise to bolster allocation
        bucketRatios[0] = 0.2;
        bucketRatios[1] = 0.4;
        bucketRatios[2] = 0.3;
        bucketRatios[3] = 0.1;

        double potentialAllocTotal = sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);
        assertEquals(550.0, potentialAllocTotal);

        allocTotal = sample.allocateBucketCounts(allocCounts, 0.5);
        assertEquals(allocTotal, potentialAllocTotal, 0.01);
        assertEquals(520.0, sample.getAllocatedCount());
        assertEquals(0.52, sample.getAllocPercent());
        assertEquals(30.0, sample.getAllocNoise());

        // finally test with ratio ranges
        for (int i = 0; i < bucketCount; ++i)
        {
            ratioRanges[i] = 0.05;
        }

        sample.clearAllocations(true);

        potentialAllocTotal = sample.getPotentialCounts(bucketRatios, requiredBuckets, ratioRanges, allocCounts);

        allocTotal = sample.allocateBucketCounts(allocCounts, 0.5);
        assertEquals(allocTotal, potentialAllocTotal, 0.01);
    }

}
