package com.hartwig.hmftools.sage.coverage;

import static com.hartwig.hmftools.sage.coverage.GeneDepth.MAX_DEPTH_BUCKET;
import static com.hartwig.hmftools.sage.coverage.GeneDepth.populateCoverageBuckets;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GeneCoverageTest
{
    @Test
    public void testMissedVariantLikelihood()
    {
        int[] baseCoverage = new int[37];
        baseCoverage[15] = 100;
        baseCoverage[20] = 100;
        baseCoverage[31] = 800;

        populateCoverageBuckets();

        assertEquals(0.002302649, GeneDepth.missedVariantLikelihood(baseCoverage), 1e-9);
    }

    @Test
    public void testMissedVariantLikelihoodNoCoverage()
    {
        int[] baseCoverage = new int[37];
        baseCoverage[0] = 1;

        populateCoverageBuckets();

        assertEquals(1, GeneDepth.missedVariantLikelihood(baseCoverage), 1e-9);
    }

    @Test
    public void testBucket()
    {
        populateCoverageBuckets();

        for(int depth = 0; depth < 30; depth++)
        {
            assertEquals(depth, GeneDepth.bucket(depth));
        }

        assertBucketAndDepth(30, 39, 30, 35);
        assertBucketAndDepth(40, 49, 31, 45);
        assertBucketAndDepth(50, 59, 32, 55);
        assertBucketAndDepth(60, 69, 33, 65);
        assertBucketAndDepth(70, 79, 34, 75);
        assertBucketAndDepth(80, 89, 35, 85);
        assertBucketAndDepth(90, 99, 36, 95);
        assertBucketAndDepth(100, 149, 37, 125);
        assertBucketAndDepth(150, 200, 38, 175);
        assertBucketAndDepth(500, 600, 45, 550);
        assertBucketAndDepth(1900, 2000, 59, 1950);
        assertBucketAndDepth(9000, 10000, 67, 9500);

        assertEquals(68, GeneDepth.bucket(10000));
        assertEquals(68, GeneDepth.bucket(MAX_DEPTH_BUCKET));
        assertEquals(10000, GeneDepth.depth(68));
        assertEquals(MAX_DEPTH_BUCKET, GeneDepth.depth(68));
    }

    private static void assertBucketAndDepth(int minDepth, int maxDepth, int expectedBucket, int expectedDepth)
    {
        for(int depth = minDepth; depth < maxDepth; depth++)
        {
            assertEquals(expectedBucket, GeneDepth.bucket(depth));
        }

        assertEquals(expectedDepth, GeneDepth.depth(expectedBucket));
    }

}
