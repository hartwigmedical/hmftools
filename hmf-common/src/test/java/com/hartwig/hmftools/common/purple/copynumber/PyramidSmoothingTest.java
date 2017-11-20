package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;

import org.junit.Test;

public class PyramidSmoothingTest {

    private final static String CHROMOSOME = "1";

    @Test
    public void testValidGermlineAmplication() {
        double neighbour  = 2.1;
        double normalRatio = 1.2;
        double expectedOKTumorCopyNumber = Math.ceil(2 * normalRatio) * neighbour;

        assertTrue(PyramidSmoothing.isValidGermlineAmplification(neighbour, createGermlineAmplified(expectedOKTumorCopyNumber + 0.01, normalRatio)));
        assertTrue(PyramidSmoothing.isValidGermlineAmplification(neighbour, createGermlineAmplified(expectedOKTumorCopyNumber, normalRatio)));
        assertFalse(PyramidSmoothing.isValidGermlineAmplification(neighbour, createGermlineAmplified(expectedOKTumorCopyNumber - 0.01, normalRatio)));
    }

    private FittedRegion createGermlineAmplified(double tumorCopyNumber, double observedNormalRatio) {
        return PurpleDatamodelTest.createDefaultFittedRegion(CHROMOSOME, 1, 2)
                .status(ObservedRegionStatus.GERMLINE_AMPLIFICATION)
                .tumorCopyNumber(tumorCopyNumber)
                .observedNormalRatio(observedNormalRatio)
                .build();
    }

}
