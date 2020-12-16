package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class LimsCohortModelFactoryTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims/cohort").getPath();

    @Test
    public void canReadCohortConfigFile() throws IOException {
        LimsCohortModel cohortModel = LimsCohortModelFactory.read(LIMS_DIRECTORY + File.separator + "cohort.tsv");

        assertNotNull(cohortModel);
        assertEquals(1, cohortModel.limsCohortMap().size());

        LimsCohortConfig cohortConfigData = cohortModel.limsCohortMap().get("CPCT");
        assertEquals("CPCT", cohortConfigData.cohortId());
        assertTrue(cohortConfigData.hospitalCenterId());
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
        assertFalse(cohortConfigData.requireAdditionalInformationForSidePanel());
    }
}