package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SampleReportFactoryTest {

    @Test
    public void canConvertHospitalPathologySampleId() {
        String wideSampleId = "WIDE020000001T";
        String cpctSampleId = "CPCT020000001T";
        String coreSampleId = "CORE020000001T";

        String correctIdT = "T20-72346";
        String correctIdC = "C18-00124";
        String wrongId = "BGr-12111";

        LimsCohortConfigData cohortConfigWIDE =
                buildTestCohortModel("WIDE", true, true, true, true, true, false, true, true, false, false, false, true);
        assertEquals(correctIdT,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT,
                        wideSampleId,
                        cohortConfigWIDE));
        assertEquals(correctIdC,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC,
                        wideSampleId,
                        cohortConfigWIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId,
                wideSampleId,
                cohortConfigWIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Lims.NOT_AVAILABLE_STRING,
                wideSampleId,
                cohortConfigWIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY,
                wideSampleId,
                cohortConfigWIDE));

        LimsCohortConfigData cohortConfigCORE =
                buildTestCohortModel("CORE", true, true, false, true, true, true, true, false, true, true, true, true);
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY,
                coreSampleId,
                cohortConfigCORE));
        assertEquals(correctIdT,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT,
                        coreSampleId,
                        cohortConfigCORE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId,
                coreSampleId,
                cohortConfigCORE));

        LimsCohortConfigData cohortConfigCPCT =
                buildTestCohortModel("CPCT", true, false, false, false, false, false, false, true, false, false, false, false);
        assertNull(correctIdT,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT,
                        cpctSampleId,
                        cohortConfigCPCT));
        assertNull(correctIdC,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC,
                        cpctSampleId,
                        cohortConfigCPCT));
    }

    @Test
    public void canCheckHospitalPatientId() {
        String coreSampleId = "CORE020000001T";
        String wideSampleId = "WIDE020000001T";

        String hospitalIdNA = Lims.NOT_AVAILABLE_STRING;
        String hospitalIDEmpty = Strings.EMPTY;
        String hospitalId = "1234";

        LimsCohortConfigData cohortConfigCORE =
                buildTestCohortModel("CORE", true, true, false, true, true, true, true, false, true, true, true, true);
        assertEquals(hospitalIdNA,
                SampleReportFactory.checkHospitalPatientId(hospitalIdNA,
                        coreSampleId,
                        cohortConfigCORE));
        assertEquals(hospitalIDEmpty,
                SampleReportFactory.checkHospitalPatientId(hospitalIDEmpty,
                        coreSampleId,
                        cohortConfigCORE));
        assertEquals(hospitalId,
                SampleReportFactory.checkHospitalPatientId(hospitalId,
                        coreSampleId,
                        cohortConfigCORE));

        LimsCohortConfigData cohortConfigWIDE =
                buildTestCohortModel("WIDE", true, true, true, true, true, false, true, true, false, false, false, true);
        assertEquals(hospitalIdNA,
                SampleReportFactory.checkHospitalPatientId(hospitalIdNA,
                        wideSampleId,
                        cohortConfigWIDE));
    }

    @NotNull
    private static LimsCohortConfigData buildTestCohortModel(@NotNull String cohortString, boolean hospitalIdCentra, boolean Report_germline,
            boolean Report_germline_flag, boolean Report_conclusion, boolean Report_viral, boolean Require_hospital_ID,
            boolean Require_hospital_PA_ID, boolean personsStudy, boolean personsrequester, boolean outputFile, boolean submission,
            boolean sidePanelInfo) {
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