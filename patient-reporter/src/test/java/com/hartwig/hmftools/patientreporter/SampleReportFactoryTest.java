package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

import org.apache.logging.log4j.util.Strings;
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

        LimsCohortConfig cohortConfigWIDE =
                PatientReportUtils.createCohortConfig("WIDE", true, true, true, true, true, false, true, true, false, true);
        assertEquals(correctIdT, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, wideSampleId, cohortConfigWIDE));
        assertEquals(correctIdC, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC, wideSampleId, cohortConfigWIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId, wideSampleId, cohortConfigWIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Lims.NOT_AVAILABLE_STRING, wideSampleId, cohortConfigWIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY, wideSampleId, cohortConfigWIDE));

        LimsCohortConfig cohortConfigCORE =
                PatientReportUtils.createCohortConfig("CORE", true, true, false, true, true, true, true, false, true, true);
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY, coreSampleId, cohortConfigCORE));
        assertEquals(correctIdT, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, coreSampleId, cohortConfigCORE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId, coreSampleId, cohortConfigCORE));

        LimsCohortConfig cohortConfigCPCT = PatientReportUtils.createCohortConfig("CPCT",
                true,
                false,
                false,
                false,
                false,
                false,
                false,
                true,
                false,
                false);
        assertNull(correctIdT, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, cpctSampleId, cohortConfigCPCT));
        assertNull(correctIdC, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC, cpctSampleId, cohortConfigCPCT));
    }

    @Test
    public void canCheckHospitalPatientId() {
        String coreSampleId = "CORE020000001T";
        String wideSampleId = "WIDE020000001T";

        String hospitalIdNA = Lims.NOT_AVAILABLE_STRING;
        String hospitalIDEmpty = Strings.EMPTY;
        String hospitalId = "1234";

        LimsCohortConfig cohortConfigCORE =
                PatientReportUtils.createCohortConfig("CORE", true, true, false, true, true, true, true, false, true, true);
        assertEquals(hospitalIdNA, SampleReportFactory.checkHospitalPatientId(hospitalIdNA, coreSampleId, cohortConfigCORE));
        assertEquals(hospitalIDEmpty, SampleReportFactory.checkHospitalPatientId(hospitalIDEmpty, coreSampleId, cohortConfigCORE));
        assertEquals(hospitalId, SampleReportFactory.checkHospitalPatientId(hospitalId, coreSampleId, cohortConfigCORE));

        LimsCohortConfig cohortConfigWIDE =
                PatientReportUtils.createCohortConfig("WIDE", true, true, true, true, true, false, true, true, false, true);
        assertEquals(hospitalIdNA, SampleReportFactory.checkHospitalPatientId(hospitalIdNA, wideSampleId, cohortConfigWIDE));
    }

}