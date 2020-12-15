package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.*;

import com.hartwig.hmftools.common.lims.LimsTestUtil;

import org.junit.Test;

public class LimsCohortModelTest {

    @Test
    public void canExtractLimsCohortModel() {
        LimsCohortConfigData cohortConfigData =
                LimsTestUtil.buildTestCohortModel("DRUP", true, false, false, false, false, false, false, true, false, false, false, false);
        assertEquals("DRUP", cohortConfigData.cohortId());
        assertTrue(cohortConfigData.hospitalCentraId());
        assertFalse(cohortConfigData.reportGermline());
        assertFalse(cohortConfigData.reportGermlineFlag());
        assertFalse(cohortConfigData.reportConclusion());
        assertFalse(cohortConfigData.reportViral());
        assertFalse(cohortConfigData.requireHospitalId());
        assertFalse(cohortConfigData.requireHospitalPAId());
        assertTrue(cohortConfigData.requireHospitalPersonsStudy());
        assertFalse(cohortConfigData.requireHospitalPersonsRequester());
        assertFalse(cohortConfigData.requirePatientIdForPdfName());
        assertFalse(cohortConfigData.requireSubmissionInformation());
        assertFalse(cohortConfigData.requireAdditionalInfromationForSidePanel());
    }
}