package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class LimsCohortModelFactoryTest {

    private static final String COHORT_MODEL_TSV = Resources.getResource("lims/cohort_config.tsv").getPath();

    @Test
    public void canReadCohortConfigFile() throws IOException {
        LimsCohortModel cohortModel = LimsCohortModelFactory.read(COHORT_MODEL_TSV);

        assertEquals(1, cohortModel.limsCohortMap().size());

        LimsCohortConfig cohortConfigData = cohortModel.limsCohortMap().get("TEST");
        assertTrue(cohortConfigData.sampleContainsHospitalCenterId());
        assertFalse(cohortConfigData.reportGermline());
        assertFalse(cohortConfigData.reportGermlineFlag());
        assertFalse(cohortConfigData.reportConclusion());
        assertFalse(cohortConfigData.reportViral());
        assertFalse(cohortConfigData.reportPeach());
        assertFalse(cohortConfigData.requireHospitalId());
        assertFalse(cohortConfigData.requireHospitalPAId());
        assertTrue(cohortConfigData.requireHospitalPersonsStudy());
        assertFalse(cohortConfigData.requireHospitalPersonsRequester());
        assertFalse(cohortConfigData.requireAdditionalInformationForSidePanel());
    }
}