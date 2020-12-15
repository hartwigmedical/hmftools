package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;

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

        assertEquals(correctIdT,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT,
                        wideSampleId,
                        buildTestCohortModelPAID("WIDE", true).queryCohortData("WIDE", wideSampleId)));
        assertEquals(correctIdC,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC,
                        wideSampleId,
                        buildTestCohortModelPAID("WIDE", true).queryCohortData("WIDE", wideSampleId)));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId,
                wideSampleId,
                buildTestCohortModelPAID("WIDE", true).queryCohortData("WIDE", wideSampleId)));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Lims.NOT_AVAILABLE_STRING,
                wideSampleId,
                buildTestCohortModelPAID("WIDE", true).queryCohortData("WIDE", wideSampleId)));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY,
                wideSampleId,
                buildTestCohortModelPAID("WIDE", true).queryCohortData("WIDE", wideSampleId)));

        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY,
                coreSampleId,
                buildTestCohortModelPAID("CORE", true).queryCohortData("CORE", coreSampleId)));
        assertEquals(correctIdT,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT,
                        coreSampleId,
                        buildTestCohortModelPAID("CORE", true).queryCohortData("CORE", coreSampleId)));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId,
                coreSampleId,
                buildTestCohortModelPAID("CORE", true).queryCohortData("CORE", coreSampleId)));

        assertNull(correctIdT,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT,
                        cpctSampleId,
                        buildTestCohortModelPAID("CPCT", false).queryCohortData("CPCT", cpctSampleId)));
        assertNull(correctIdC,
                SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC,
                        cpctSampleId,
                        buildTestCohortModelPAID("CPCT", false).queryCohortData("CPCT", cpctSampleId)));
    }

    @Test
    public void canCheckHospitalPatientId() {
        String coreSampleId = "CORE020000001T";
        String wideSampleId = "WIDE020000001T";

        String hospitalIdNA = Lims.NOT_AVAILABLE_STRING;
        String hospitalIDEmpty = Strings.EMPTY;
        String hospitalId = "1234";

        assertEquals(hospitalIdNA,
                SampleReportFactory.checkHospitalPatientId(hospitalIdNA,
                        coreSampleId,
                        buildTestCohortModelPatientId("CORE", true).queryCohortData("CORE", coreSampleId)));
        assertEquals(hospitalIDEmpty,
                SampleReportFactory.checkHospitalPatientId(hospitalIDEmpty,
                        coreSampleId,
                        buildTestCohortModelPatientId("CORE", true).queryCohortData("CORE", coreSampleId)));
        assertEquals(hospitalId,
                SampleReportFactory.checkHospitalPatientId(hospitalId,
                        coreSampleId,
                        buildTestCohortModelPatientId("CORE", true).queryCohortData("CORE", coreSampleId)));
        assertEquals(hospitalIdNA,
                SampleReportFactory.checkHospitalPatientId(hospitalIdNA,
                        wideSampleId,
                        buildTestCohortModelPatientId("WIDE", true).queryCohortData("WIDE", wideSampleId)));
    }

    @NotNull
    private static LimsCohortModel buildTestCohortModelPAID(@NotNull String cohortString, boolean requirePAId) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(requirePAId)
                .hospitalPersonsStudy(true)
                .hospitalPersonsRequester(false)
                .outputFile(false)
                .submission(false)
                .sidePanelInfo(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }

    @NotNull
    private static LimsCohortModel buildTestCohortModelPatientId(@NotNull String cohortString, boolean requirePatientId) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(requirePatientId)
                .requireHospitalPAId(false)
                .hospitalPersonsStudy(true)
                .hospitalPersonsRequester(false)
                .outputFile(false)
                .submission(false)
                .sidePanelInfo(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }
}