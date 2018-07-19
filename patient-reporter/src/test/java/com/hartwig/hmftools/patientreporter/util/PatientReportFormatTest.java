package com.hartwig.hmftools.patientreporter.util;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;

import org.junit.Test;

public class PatientReportFormatTest {

    @Test
    public void removesCopiesWhenNoTumor() {
        assertTrue(PatientReportFormat.correctValueForFitStatus(FittedPurityStatus.NORMAL, "1").equals("1"));
        assertFalse(PatientReportFormat.correctValueForFitStatus(FittedPurityStatus.NO_TUMOR, "1").equals("1"));
    }
}