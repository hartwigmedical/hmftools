package com.hartwig.hmftools.cobalt.consolidation;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class ResultsConsolidatorTest
{
    @Test
    public void testCalcConsolidationCount()
    {
        assertEquals(1, ResultsConsolidator.calcConsolidationCount(16.0));
        assertEquals(1, ResultsConsolidator.calcConsolidationCount(8.1));

        // starting from 8.0 reads depth we consolidate
        assertEquals(10, ResultsConsolidator.calcConsolidationCount(8.0));
        assertEquals(30, ResultsConsolidator.calcConsolidationCount(3.0));
        assertEquals(100, ResultsConsolidator.calcConsolidationCount(0.8));
        assertEquals(300, ResultsConsolidator.calcConsolidationCount(0.3));
        assertEquals(500, ResultsConsolidator.calcConsolidationCount(0.15));
        assertEquals(1000, ResultsConsolidator.calcConsolidationCount(0.08));
        assertEquals(1000, ResultsConsolidator.calcConsolidationCount(0.008));
        assertEquals(1000, ResultsConsolidator.calcConsolidationCount(0.00015));
    }

    @Test
    public void testCalcConsolidateBoundaries()
    {
        List<Integer> windowPositions = new ArrayList<>();
        windowPositions.add(1_001);
        windowPositions.add(2_001);
        windowPositions.add(5_001);
        windowPositions.add(6_001);
        windowPositions.add(10_001);

        // test we do not span over centromere
        windowPositions.add(3_020_001);
        windowPositions.add(3_021_001);
        windowPositions.add(3_024_001);
        windowPositions.add(3_029_001);

        List<LowCovBucket> buckets = ResultsConsolidator.consolidateIntoBuckets(windowPositions, 4);

        assertEquals(3, buckets.size());

        // check each one
        assertEquals(1_001, buckets.get(0).StartPosition);
        assertEquals(8_001, buckets.get(0).EndPosition);
        assertEquals(9_001, buckets.get(1).StartPosition);
        assertEquals(11_001, buckets.get(1).EndPosition);
        assertEquals(3_020_001, buckets.get(2).StartPosition);
        assertEquals(3_030_001, buckets.get(2).EndPosition);

        windowPositions.add(3_034_001);

        buckets = ResultsConsolidator.consolidateIntoBuckets(windowPositions, 4);

        assertEquals(4, buckets.size());

        // check each one
        assertEquals(1_001, buckets.get(0).StartPosition);
        assertEquals(8_001, buckets.get(0).EndPosition);
        assertEquals(9_001, buckets.get(1).StartPosition);
        assertEquals(11_001, buckets.get(1).EndPosition);
        assertEquals(3_020_001, buckets.get(2).StartPosition);
        assertEquals(3_031_001, buckets.get(2).EndPosition);
        assertEquals(3_032_001, buckets.get(3).StartPosition);
        assertEquals(3_035_001, buckets.get(3).EndPosition);
    }
}
