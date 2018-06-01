package com.hartwig.hmftools.patientreporter.copynumber;

import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.ABS_LOSS;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.REL_GAIN;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.includeInReport;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ReportableCopyNumbersTest {

    @Test
    public void testLoss() {
        assertTrue(includeInReport(0, ABS_LOSS));
        assertTrue(includeInReport(1, ABS_LOSS));
        assertTrue(includeInReport(2, ABS_LOSS));
        assertTrue(includeInReport(2, ABS_LOSS - 0.1));
        assertFalse(includeInReport(2, ABS_LOSS + 0.1));
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
