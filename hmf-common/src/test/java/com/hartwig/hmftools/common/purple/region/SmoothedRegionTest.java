package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SmoothedRegionTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void stuff() {
        assertEquals(1.25, SmoothedRegions.allowedCopyNumberDeviation(0), EPSILON);
        assertEquals(0.75, SmoothedRegions.allowedCopyNumberDeviation(5), EPSILON);
        assertEquals(0.25, SmoothedRegions.allowedCopyNumberDeviation(10), EPSILON);
        assertEquals(0.25, SmoothedRegions.allowedCopyNumberDeviation(50), EPSILON);
    }

}
