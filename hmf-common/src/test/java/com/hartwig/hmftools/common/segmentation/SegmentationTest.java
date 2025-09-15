package com.hartwig.hmftools.common.segmentation;

import org.junit.Test;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class SegmentationTest extends SegmentationTestBase
{
    @Test
    public void completeSeries() {
        assertArrayEquals(d(1.0, 1.0), segmentation(d(1.0, 1.0)).completeSeries(), 1e-10);
        assertArrayEquals(d(1.0, 1.0, 2.0, 2.1), segmentation(d(1.0, 1.0), d(2.0, 2.1)).completeSeries(), 1e-10);
        assertArrayEquals(d(10.0, 1.0, 2.0, 2.1), segmentation(d(10.0), d(1.0), d(2.0, 2.1)).completeSeries(), 1e-10);
    }

    @Test
    public void cost() {
        assertEquals(300.0, segmentation(d(2.0), d(3.0), d(1.0), d(12.0), d(13.0), d(11.0)).cost(50.0), 1e-10);
        assertEquals(52.0, segmentation(d(2.0, 3.0, 1.0)).cost(50.0), 1e-10);
        assertEquals(104.0, segmentation(d(2.0, 3.0, 1.0), d(12.0, 13.0, 11.0)).cost(50.0), 1e-10);

        assertEquals(300.0, segmentation(d(2), d(3), d(1), d(12), d(13), d(11)).cost(50.0), 1e-10);
        assertEquals(250.5, segmentation(d(2, 3), d(1), d(12), d(13), d(11)).cost(50.0), 1e-10);
        assertEquals(202.0, segmentation(d(2, 3, 1), d(12), d(13), d(11)).cost(50.0), 1e-10);
        assertEquals(152.5, segmentation(d(2, 3, 1), d(12, 13), d(11)).cost(50.0), 1e-10);
        assertEquals(104.0, segmentation(d(2, 3, 1), d(12, 13, 11)).cost(50.0), 1e-10);
        assertEquals(179.0, segmentation(d(2, 3, 1, 12), d(13, 11)).cost(50.0), 1e-10);
        assertEquals(213.0, segmentation(d(2, 3), d(1, 12), d(13, 11)).cost(50.0), 1e-10);
        assertEquals(193.25, segmentation(d(2, 3), d(1, 12, 13, 11)).cost(50.0), 1e-10);
        assertEquals(227.0, segmentation(d(2, 3, 1, 12), d(13), d(11)).cost(50.0), 1e-10);
        assertEquals(204.0, segmentation(d(2, 3, 1, 12, 13, 11)).cost(50.0), 1e-10);
        assertEquals(252.0, segmentation(d(2), d(3), d(1), d(12), d(13, 11)).cost(50.0), 1e-10);
    }
}