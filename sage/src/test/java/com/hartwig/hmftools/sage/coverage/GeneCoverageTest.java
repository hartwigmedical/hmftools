package com.hartwig.hmftools.sage.coverage;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GeneCoverageTest {

    @Test
    public void testMissedVariantLikelihood() {
        int[] baseCoverage = new int[37];
        baseCoverage[15] = 100;
        baseCoverage[20] = 100;
        baseCoverage[31] = 800;

        assertEquals(0.002302649, GeneCoverage.missedVariantLikelihood(baseCoverage), 1e-9);
    }

    @Test
    public void testMissedVariantLikelihoodNoCoverage() {
        int[] baseCoverage = new int[37];
        baseCoverage[0] = 1;

        assertEquals(1, GeneCoverage.missedVariantLikelihood(baseCoverage), 1e-9);
    }

    @Test
    public void testBucket() {
        for (int depth = 0; depth < 30; depth++) {
            assertEquals(depth, GeneCoverage.bucket(depth));
        }

        assertBucketAndDepth(30, 39, 30, 35);
        assertBucketAndDepth(40, 49, 31, 45);
        assertBucketAndDepth(50, 59, 32, 55);
        assertBucketAndDepth(60, 69, 33, 65);
        assertBucketAndDepth(70, 79, 34, 75);
        assertBucketAndDepth(80, 89, 35, 85);
        assertBucketAndDepth(90, 99, 36, 95);
        assertBucketAndDepth(100, 1000, 37, 100);
    }

    private static void assertBucketAndDepth(int minDepth, int maxDepth, int expectedBucket, int expectedDepth) {
        for (int depth = minDepth; depth < maxDepth; depth++) {
            assertEquals(expectedBucket, GeneCoverage.bucket(depth));
        }

        assertEquals(expectedDepth, GeneCoverage.depth(expectedBucket));
    }

}
