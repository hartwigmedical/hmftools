package com.hartwig.hmftools.patientreporter.report.util;

import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;

import org.junit.Test;

public class PatientReportFormatTest {

    @Test
    public void removesCopiesWhenNoTumor() {
        assertTrue(PatientReportFormat.correctValueForFitStatus(FittedPurityStatus.NORMAL, "1").equals("1"));

        assertTrue(PatientReportFormat.correctValueForFitStatus(FittedPurityStatus.NO_TUMOR, "1").equals("N/A"));
    }
}