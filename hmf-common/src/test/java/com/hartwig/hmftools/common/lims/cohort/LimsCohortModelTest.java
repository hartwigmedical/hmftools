package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.*;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LimsCohortModelTest {

    @Test
    public void canExtractLimsCohortModel() {
        LimsCohortConfigData cohortConfigData = buildTestCohortModel("DRUP").queryCohortData("DRUP", "DRUP02010001T");
        assertEquals("DRUP", cohortConfigData.cohortId());
        assertTrue(cohortConfigData.hospitalCentraId());
        assertFalse(cohortConfigData.reportGermline());
        assertFalse(cohortConfigData.reportGermlineFlag());
        assertFalse(cohortConfigData.reportConclusion());
        assertFalse(cohortConfigData.reportViral());
        assertFalse(cohortConfigData.requireHospitalId());
        assertFalse(cohortConfigData.requireHospitalPAId());
        assertFalse(cohortConfigData.requireHospitalPersonsStudy());
        assertFalse(cohortConfigData.requireHospitalPersonsRequester());
        assertFalse(cohortConfigData.requirePatientIdForPdfName());
        assertFalse(cohortConfigData.requireSubmissionInformation());
        assertFalse(cohortConfigData.requireAdditionalInfromationForSidePanel());
    }

    @NotNull
    private static LimsCohortModel buildTestCohortModel(@NotNull String cohortString) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalCentraId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .requireHospitalPersonsStudy(false)
                .requireHospitalPersonsRequester(false)
                .requirePatientIdForPdfName(false)
                .requireSubmissionInformation(false)
                .requireAdditionalInfromationForSidePanel(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }

}