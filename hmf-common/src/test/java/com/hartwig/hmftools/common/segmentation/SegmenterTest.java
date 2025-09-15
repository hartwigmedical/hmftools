package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import org.junit.Test;

public class SegmenterTest extends SegmentationTestBase
{
    @Test
    public void compareAlgorithmWithExhaustiveSearch()
    {
        compareEvaluationMethods(d(1, 1, 2, 2), 50.0);
        compareEvaluationMethods(d(1, 1, 2, 2), 10.0);
        compareEvaluationMethods(d(-1, 1, -2, 2), 10.0);
        compareEvaluationMethods(d(2, 3, 1, 12, 13, 11), 50.0);
        compareEvaluationMethods(d(5, 4, 1, 1), 0.0);
        compareEvaluationMethods(d(5, 4, 1, 1), 50.0);
        compareEvaluationMethods(d(1, 1, -5, -6, 1, 2), 50.0);
        compareEvaluationMethods(d(55, 55, -276, -275, 55, 56), 50.0);
    }

    private void compareEvaluationMethods(double[] data, double gamma)
    {
        Segmenter segmenter = new Segmenter(data, gamma, true);
        double segmentCost = segmenter.segmentPenalty;
        double leastCostByFastSearch = segmenter.leastCostSegmentation.cost(segmentCost);
        ExhaustiveSearchSegmenter exhaustiveSearchSegmenter = new ExhaustiveSearchSegmenter(data, gamma, true);
        double leastCostByExhaustiveSearch = exhaustiveSearchSegmenter.cheapestSegmentationByExhaustiveSearch().cost(segmentCost);
        assertEquals(leastCostByFastSearch, leastCostByExhaustiveSearch, 0.001);
    }

    @Test
    public void segmentBy()
    {
        assertEquals(
                segmentation(d(2), d(3), d(1), d(12), d(13), d(11)),
                new Segmenter(d(2, 3, 1, 12, 13, 11)).segmentBy(java.util.Arrays.asList(0, 1, 2, 3, 4, 5))
        );
        assertEquals(
                segmentation(d(2), d(3), d(1), d(12), d(13, 11)),
                new Segmenter(d(2, 3, 1, 12, 13, 11)).segmentBy(java.util.Arrays.asList(0, 1, 2, 3, 5))
        );
        assertEquals(
                segmentation(d(2), d(3, 1, 12), d(13, 11)),
                new Segmenter(d(2, 3, 1, 12, 13, 11)).segmentBy(java.util.Arrays.asList(0, 3, 5))
        );
        assertEquals(
                segmentation(d(2), d(3, 1, 12, 13, 11)),
                new Segmenter(d(2, 3, 1, 12, 13, 11)).segmentBy(java.util.Arrays.asList(0, 5))
        );
        assertEquals(
                segmentation(d(2, 3, 1, 12, 13, 11)),
                new Segmenter(d(2, 3, 1, 12, 13, 11)).segmentBy(List.of(5))
        );
    }

    @Test
    public void pcf()
    {
        assertEquals(
                new PiecewiseConstantFit(new int[] { 1 }, new int[] { 0 }, d(1.0)),
                new Segmenter(d(1.0), 50.0, false).pcf()
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 2 }, new int[] { 0 }, d(1.0)),
                new Segmenter(d(1.0, 1.0), 50.0, false).pcf()
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 3 }, new int[] { 0 }, d(1.0)),
                new Segmenter(d(1.0, 1.0, 1.0), 50.0, false).pcf()
        );

        assertEquals(
                new PiecewiseConstantFit(new int[] { 3, 3 }, new int[] { 0, 3 }, d(2.0, 12.0)),
                new Segmenter(d(2, 3, 1, 12, 13, 11), 50.0, false).pcf()
        );

        assertEquals(
                new PiecewiseConstantFit(new int[] { 3, 3, 3 }, new int[] { 0, 3, 6 }, d(2.0, 12.0, 2.0)),
                new Segmenter(d(2, 3, 1, 12, 13, 11, 2, 3, 1), 50.0, false).pcf()
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 3, 3, 3 }, new int[] { 0, 3, 6 }, d(2.0, 12.0, 22.0)),
                new Segmenter(d(2, 3, 1, 12, 13, 11, 22, 23, 21), 50.0, false).pcf()
        );
    }

    // Tests for ExhaustiveSearchSegmenter so that we know it can be trusted in tests.
    @Test
    public void exhaustiveSearchSegmenterTest()
    {
        Set<Segmentation> segmentations1 = new ExhaustiveSearchSegmenter(d(1)).allPossibleSegmentations();
        assertEquals(1, segmentations1.size());
        assertTrue(segmentations1.contains(segmentation(d(1))));

        Set<Segmentation> segmentations2 = new ExhaustiveSearchSegmenter(d(1, 1)).allPossibleSegmentations();
        assertEquals(2, segmentations2.size());
        assertTrue(segmentations2.contains(segmentation(d(1), d(1))));
        assertTrue(segmentations2.contains(segmentation(d(1, 1))));

        Set<Segmentation> segmentations3 = new ExhaustiveSearchSegmenter(d(1, 2, 3)).allPossibleSegmentations();
        assertEquals(4, segmentations3.size());
        assertTrue(segmentations3.contains(segmentation(d(1), d(2), d(3))));
        assertTrue(segmentations3.contains(segmentation(d(1, 2), d(3))));
        assertTrue(segmentations3.contains(segmentation(d(1), d(2, 3))));
        assertTrue(segmentations3.contains(segmentation(d(1, 2, 3))));

        Set<Segmentation> segmentations4 = new ExhaustiveSearchSegmenter(d(1, 2, 3, 1)).allPossibleSegmentations();
        assertEquals(8, segmentations4.size());
        assertTrue(segmentations4.contains(segmentation(d(1), d(2), d(3), d(1))));
        assertTrue(segmentations4.contains(segmentation(d(1, 2), d(3), d(1))));
        assertTrue(segmentations4.contains(segmentation(d(1), d(2, 3), d(1))));
        assertTrue(segmentations4.contains(segmentation(d(1), d(2), d(3, 1))));
        assertTrue(segmentations4.contains(segmentation(d(1, 2), d(3, 1))));
        assertTrue(segmentations4.contains(segmentation(d(1), d(2, 3, 1))));
        assertTrue(segmentations4.contains(segmentation(d(1, 2, 3), d(1))));
        assertTrue(segmentations4.contains(segmentation(d(1, 2, 3, 1))));
    }

    @Test
    public void exhaustiveSearchCheapestTest()
    {
        assertEquals(
                new PiecewiseConstantFit(new int[] { 1 }, new int[] { 0 }, d(1.0)),
                new ExhaustiveSearchSegmenter(d(1.0), 50.0, false).cheapestSegmentationByExhaustiveSearch().pcf()
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 2 }, new int[] { 0 }, d(1.0)),
                new ExhaustiveSearchSegmenter(d(1.0, 1.0), 50.0, false).cheapestSegmentationByExhaustiveSearch().pcf()
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 3 }, new int[] { 0 }, d(1.0)),
                new ExhaustiveSearchSegmenter(d(1.0, 1.0, 1.0), 50.0, false).cheapestSegmentationByExhaustiveSearch().pcf()
        );

        assertEquals(
                new PiecewiseConstantFit(new int[] { 3, 3 }, new int[] { 0, 3 }, d(2.0, 12.0)),
                new ExhaustiveSearchSegmenter(d(2, 3, 1, 12, 13, 11), 50.0, false).cheapestSegmentationByExhaustiveSearch().pcf()
        );

        assertEquals(
                new PiecewiseConstantFit(new int[] { 3, 3, 3 }, new int[] { 0, 3, 6 }, d(2.0, 12.0, 2.0)),
                new ExhaustiveSearchSegmenter(d(2, 3, 1, 12, 13, 11, 2, 3, 1), 50.0, false).cheapestSegmentationByExhaustiveSearch().pcf()
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 3, 3, 3 }, new int[] { 0, 3, 6 }, d(2.0, 12.0, 22.0)),
                new ExhaustiveSearchSegmenter(d(2, 3, 1, 12, 13, 11, 22, 23, 21), 50.0, false).cheapestSegmentationByExhaustiveSearch()
                        .pcf()
        );
    }
}
