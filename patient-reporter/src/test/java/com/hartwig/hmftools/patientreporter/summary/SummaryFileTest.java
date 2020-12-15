package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SummaryFileTest {

    private static final String SAMPLE_SUMMARY_TSV = Resources.getResource("sample_summary/sample_summary.tsv").getPath();

    @Test
    public void summaryFromCSVWithNewLines() throws IOException {
        SummaryModel summaryModel = SummaryFile.buildFromTsv(SAMPLE_SUMMARY_TSV);
        assertEquals(1, summaryModel.summaryCount());

        LimsCohortConfigData cohortConfig =
                buildTestCohortModel("CORE", true, true, false, true, true, true, true, false, true, true, true, true);

        String summary = summaryModel.findSummaryForSample("sample", cohortConfig);

        assertEquals(3, summary.split("\n").length);
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