package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SummaryModelTest {

    @Test
    public void sampleArePresentInSummaryModel() {
        Map<String, String> summaryToSampleMap = Maps.newHashMap();
        summaryToSampleMap.put("sample", "this is a test summary");
        SummaryModel summaryModel = new SummaryModel(summaryToSampleMap);

        assertTrue(summaryModel.samplePresentInSummaries("sample"));
        assertFalse(summaryModel.samplePresentInSummaries("sample2"));
    }

    @Test
    public void canExtractSummaryOfSample() {
        Map<String, String> summaryToSampleMap = Maps.newHashMap();
        summaryToSampleMap.put("sample", "this is a test summary");
        SummaryModel summaryModel = new SummaryModel(summaryToSampleMap);

        assertEquals("this is a test summary", summaryModel.findSummaryForSample("sample", buildTestCohortModel("WIDE", true).queryCohortData("WIDE", "WIDE020000001T")));
        assertNotEquals("this is a test summary", summaryModel.findSummaryForSample("sample2", buildTestCohortModel("WIDE", true).queryCohortData("WIDE", "WIDE020000001T")));
    }

    @NotNull
    private static LimsCohortModel buildTestCohortModel(@NotNull String cohortString, boolean requireConclusion) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalCentraId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(requireConclusion)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .requireHospitalPersonsStudy(true)
                .requireHospitalPersonsRequester(false)
                .requirePatientIdForPdfName(false)
                .requireSubmissionInformation(false)
                .requireAdditionalInfromationForSidePanel(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }
}