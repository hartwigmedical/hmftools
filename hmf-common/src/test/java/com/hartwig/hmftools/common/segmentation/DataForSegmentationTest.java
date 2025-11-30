package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class DataForSegmentationTest
{
    @Test
    public void countTest()
    {
        double[] segmentationValues = new double[]{1.0, 2.0, 3.0};
        double[] rawValues = new double[]{1.1, 2.1, 3.1};
        DataForSegmentation sd = new DataForSegmentation(segmentationValues, rawValues);
        assertEquals(3, sd.count());
    }

    @Test
    public void isEmptyTest()
    {
        double[] segmentationValues = new double[]{};
        double[] rawValues = new double[]{};
        DataForSegmentation sd = new DataForSegmentation(segmentationValues, rawValues);
        assertTrue(sd.isEmpty());
        segmentationValues = new double[]{1.0, 2.0, 3.0};
        rawValues = new double[]{1.1, 2.1, 3.1};
        sd = new DataForSegmentation(segmentationValues, rawValues);
        assertFalse(sd.isEmpty());
    }
}
