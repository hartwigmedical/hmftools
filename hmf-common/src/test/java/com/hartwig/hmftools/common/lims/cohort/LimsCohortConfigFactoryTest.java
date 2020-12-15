package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class LimsCohortConfigFactoryTest {
    private static final String LIMS_DIRECTORY = Resources.getResource("lims/cohort").getPath();

    @Test
    public void canReadCohortConfigFile() throws IOException {
        LimsCohortModel cohortModel = LimsCohortConfigFactory.read(LIMS_DIRECTORY + File.separator + "cohort.tsv");

        assertNotNull(cohortModel);
        assertEquals(1, cohortModel.limsCohortMap().size());

        LimsCohortConfigData cohortConfigData = cohortModel.limsCohortMap().get("CPCT");
        assertEquals("CPCT", cohortConfigData.cohortId());
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