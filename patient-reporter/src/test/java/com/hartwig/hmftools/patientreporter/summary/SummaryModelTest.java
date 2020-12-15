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

        LimsCohortConfigData cohortConfig = buildTestCohortModel("WIDE", true, true, true, true, true, false, true, true, false, false, false, true);

        assertEquals("this is a test summary",
                summaryModel.findSummaryForSample("sample", cohortConfig));
        assertNotEquals("this is a test summary",
                summaryModel.findSummaryForSample("sample2", cohortConfig));
    }

    @NotNull
    private static LimsCohortConfigData buildTestCohortModel(@NotNull String cohortString, boolean hospitalIdCentra,
            boolean Report_germline, boolean Report_germline_flag, boolean Report_conclusion, boolean Report_viral,
            boolean Require_hospital_ID, boolean Require_hospital_PA_ID, boolean personsStudy, boolean personsrequester, boolean outputFile,
            boolean submission, boolean sidePanelInfo) {
        return ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalCentraId(hospitalIdCentra)
                .reportGermline(Report_germline)
                .reportGermlineFlag(Report_germline_flag)
                .reportConclusion(Report_conclusion)
                .reportViral(Report_viral)
                .requireHospitalId(Require_hospital_ID)
                .requireHospitalPAId(Require_hospital_PA_ID)
                .requireHospitalPersonsStudy(personsStudy)
                .requireHospitalPersonsRequester(personsrequester)
                .requirePatientIdForPdfName(outputFile)
                .requireSubmissionInformation(submission)
                .requireAdditionalInfromationForSidePanel(sidePanelInfo)
                .build();
    }
}