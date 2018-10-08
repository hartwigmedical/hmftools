package com.hartwig.hmftools.patientreporter.copynumber;

import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.ABS_LOSS;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.REL_GAIN;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.includeInReport;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ReportableCopyNumbersTest {

    @Test
    public void reportLossesCorrectly() {
        assertTrue(includeInReport(2, ABS_LOSS, true, true));
        assertTrue(includeInReport(2, ABS_LOSS, false, true));
        assertFalse(includeInReport(2, ABS_LOSS, true, false));

        assertTrue(includeInReport(1, ABS_LOSS, true, true));
        assertTrue(includeInReport(2, ABS_LOSS, true, true));
        assertTrue(includeInReport(2, ABS_LOSS - 0.1, true, true));
        assertFalse(includeInReport(2, ABS_LOSS + 0.1, true, true));
    }

    @Test
    public void reportGainsCorrectly() {
        assertTrue(includeInReport(0, REL_GAIN, true, true));

        assertTrue(includeInReport(1, REL_GAIN, true, true));
        assertTrue(includeInReport(1, REL_GAIN, true, false));
        assertFalse(includeInReport(1, REL_GAIN, false, true));

        assertTrue(includeInReport(1, REL_GAIN + 0.1, true, true));
        assertFalse(includeInReport(1, REL_GAIN - 0.1, true, true));
        assertFalse(includeInReport(1.01, REL_GAIN, true, true));
    }
}
