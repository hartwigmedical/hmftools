package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.metrics.GeneCoverage.MAX_DEPTH_BUCKET;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.gene.GeneRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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

        assertEquals(0.002302649, GeneCoverage.missedVariantLikelihood(baseCoverage), 1e-9);
    }

    @Test
    public void testMissedVariantLikelihoodNoCoverage()
    {
        int[] baseCoverage = new int[37];
        baseCoverage[0] = 1;

        assertEquals(1, GeneCoverage.missedVariantLikelihood(baseCoverage), 1e-9);
    }

    @Test
    public void testBucket()
    {
        for(int depth = 0; depth < 30; depth++)
        {
            assertEquals(depth, GeneCoverage.bucket(depth));
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

        assertEquals(68, GeneCoverage.bucket(10000));
        assertEquals(68, GeneCoverage.bucket(MAX_DEPTH_BUCKET));
        assertEquals(10000, GeneCoverage.depth(68));
        assertEquals(MAX_DEPTH_BUCKET, GeneCoverage.depth(68));
    }

    private static void assertBucketAndDepth(int minDepth, int maxDepth, int expectedBucket, int expectedDepth)
    {
        for(int depth = minDepth; depth < maxDepth; depth++)
        {
            assertEquals(expectedBucket, GeneCoverage.bucket(depth));
        }

        assertEquals(expectedDepth, GeneCoverage.depth(expectedBucket));
    }

    @Test
    public void testAlignmentBefore()
    {
        ExonCoverage victim = exon("Gene", 100, 104);
        victim.processRead(30, 99);
        assertCoverage(victim, 0, 0, 0, 0, 0);
    }

    @Test
    public void testAlignmentOverlapsStart()
    {
        ExonCoverage victim = exon("Gene", 100, 104);
        victim.processRead(90, 100);
        assertCoverage(victim, 1, 0, 0, 0, 0);

        victim.processRead(91, 102);
        assertCoverage(victim, 2, 1, 1, 0, 0);
    }

    @Test
    public void testAlignmentOverlapsEnd()
    {
        ExonCoverage victim = exon("Gene", 100, 104);
        victim.processRead(103, 110);
        assertCoverage(victim, 0, 0, 0, 1, 1);

        victim.processRead(104, 110);
        assertCoverage(victim, 0, 0, 0, 1, 2);
    }

    @Test
    public void testAlignmentOverlapsBoth()
    {
        ExonCoverage victim = exon("Gene", 100, 104);
        victim.processRead(90, 110);
        assertCoverage(victim, 1, 1, 1, 1, 1);

        victim.processRead(90, 110);
        assertCoverage(victim, 2, 2, 2, 2, 2);
    }

    @Test
    public void testAlignmentWithin()
    {
        ExonCoverage victim = exon("Gene", 100, 104);
        victim.processRead(101, 103);
        assertCoverage(victim, 0, 1, 1, 1, 0);

        victim.processRead(101, 105);
        assertCoverage(victim, 0, 2, 2, 2, 1);

        victim.processRead(99, 101);
        assertCoverage(victim, 1, 3, 2, 2, 1);
    }

    @Test
    public void testAlignmentAfter()
    {
        ExonCoverage victim = exon("Gene", 100, 104);
        victim.processRead(106, 200);
        assertCoverage(victim, 0, 0, 0, 0, 0);
    }

    private void assertCoverage(ExonCoverage victim, int... values)
    {
        for(int i = 0; i < values.length; i++)
        {
            int value = values[i];
            assertEquals(value, victim.coverage()[i]);
        }
    }

    private ChrBaseRegion alignment(int start, int end)
    {
        return new ChrBaseRegion("1", start, end);
    }

    private ExonCoverage exon(String gene, int start, int end)
    {
        return new ExonCoverage(new GeneRegion(CHR_1, start, end, GENE_NAME_1, 1));
    }
}
