package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;
import static com.hartwig.hmftools.common.numeric.Doubles.lessThan;
import static com.hartwig.hmftools.common.purple.region.SmoothedRegions.allowedBAFDeviation;
import static com.hartwig.hmftools.common.purple.region.SmoothedRegions.allowedCopyNumberDeviation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SmoothedRegionTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void copyNumberAllowances() {
        assertEquals(1.3, allowedCopyNumberDeviation(0), EPSILON);
        assertEquals(0.8, allowedCopyNumberDeviation(5), EPSILON);
        assertEquals(0.3, allowedCopyNumberDeviation(10), EPSILON);
        assertEquals(0.3, allowedCopyNumberDeviation(50), EPSILON);
    }

    @Test
    public void bafAllowances() {
        assertTrue(lessThan(allowedBAFDeviation(1), 0.44));
        assertTrue(greaterThan(allowedBAFDeviation(7), 0.072));
        assertTrue(greaterThan(allowedBAFDeviation(14), 0.052));
        assertTrue(greaterThan(allowedBAFDeviation(46), 0.03));
    }
}
