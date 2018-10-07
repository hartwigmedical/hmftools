package com.hartwig.hmftools.patientreporter.copynumber;

import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.ABS_LOSS;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.REL_GAIN;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportableCopyNumbers.includeInReport;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.junit.Test;

public class ReportableCopyNumbersTest {

    @Test
    public void reportLossesCorrectly() {
        assertTrue(includeInReport(2, ABS_LOSS, DriverCategory.TSG));
        assertTrue(includeInReport(2, ABS_LOSS, null));
        assertFalse(includeInReport(2, ABS_LOSS, DriverCategory.ONCO));

        assertTrue(includeInReport(1, ABS_LOSS, DriverCategory.TSG));
        assertTrue(includeInReport(2, ABS_LOSS, DriverCategory.TSG));
        assertTrue(includeInReport(2, ABS_LOSS - 0.1, DriverCategory.TSG));
        assertFalse(includeInReport(2, ABS_LOSS + 0.1, DriverCategory.TSG));
    }

    @Test
    public void reportGainsCorrectly() {
        assertTrue(includeInReport(0, REL_GAIN, DriverCategory.ONCO));

        assertTrue(includeInReport(1, REL_GAIN, DriverCategory.ONCO));
        assertTrue(includeInReport(1, REL_GAIN, null));
        assertFalse(includeInReport(1, REL_GAIN, DriverCategory.TSG));

        assertTrue(includeInReport(1, REL_GAIN + 0.1, DriverCategory.ONCO));
        assertFalse(includeInReport(1, REL_GAIN - 0.1, DriverCategory.ONCO));
        assertFalse(includeInReport(1.01, REL_GAIN, DriverCategory.ONCO));
    }
}
