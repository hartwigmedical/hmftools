package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.LimsTestUtil;

import org.junit.Test;

public class LimsCohortModelTest {

    @Test
    public void canExtractLimsCohortModel() {
        LimsCohortConfig cohortConfigData =
                LimsTestUtil.buildTestCohortModel("DRUP", true, false, false, false, false, false, false, true, false, false, false, false);
        assertEquals("DRUP", cohortConfigData.cohortId());
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

    @Test
    public void canQueryLimsCohortModel() {
        LimsCohortConfig cohortConfigData =
                LimsTestUtil.buildTestCohortModel("DRUP", true, false, false, false, false, false, false, true, false, false, false, false);
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", cohortConfigData);
        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertNull(model.queryCohortData(null, "DRUP01"));
        assertNull(model.queryCohortData("CPCT", "DRUP01"));
        assertNull(model.queryCohortData("DRUP", "CPCT01"));
        assertEquals(cohortConfigData, model.queryCohortData("DRUP", "DRUP01"));

    }
}