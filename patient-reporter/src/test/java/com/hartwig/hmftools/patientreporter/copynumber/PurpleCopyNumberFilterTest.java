package com.hartwig.hmftools.patientreporter.copynumber;

import static com.hartwig.hmftools.patientreporter.copynumber.PurpleCopyNumberFilter.ABS_GAIN;
import static com.hartwig.hmftools.patientreporter.copynumber.PurpleCopyNumberFilter.ABS_LOSS;
import static com.hartwig.hmftools.patientreporter.copynumber.PurpleCopyNumberFilter.REL_GAIN;
import static com.hartwig.hmftools.patientreporter.copynumber.PurpleCopyNumberFilter.includeInReport;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class PurpleCopyNumberFilterTest {

    @Test
    public void testLoss() {
        assertTrue(includeInReport(0, ABS_LOSS));
        assertTrue(includeInReport(1, ABS_LOSS));
        assertTrue(includeInReport(2, ABS_LOSS));
        assertTrue(includeInReport(2, ABS_LOSS - 0.1));
        assertFalse(includeInReport(2, ABS_LOSS + 0.1));
    }

    @Test
    public void testAbsGain() {
        assertTrue(includeInReport(0, ABS_GAIN));
        assertTrue(includeInReport(1, ABS_GAIN));
        assertTrue(includeInReport(2, ABS_GAIN));
        assertTrue(includeInReport(ABS_GAIN, ABS_GAIN));

        assertTrue(includeInReport(ABS_GAIN, ABS_GAIN + 0.1));
        assertFalse(includeInReport(ABS_GAIN, ABS_GAIN - 0.1));
    }

    @Test
    public void testRelGain() {
        assertTrue(includeInReport(0, REL_GAIN));
        assertTrue(includeInReport(1, REL_GAIN));
        assertTrue(includeInReport(1, REL_GAIN + 0.1));
        assertFalse(includeInReport(1, REL_GAIN - 0.1));
        assertFalse(includeInReport(1.01, REL_GAIN));
    }
}
