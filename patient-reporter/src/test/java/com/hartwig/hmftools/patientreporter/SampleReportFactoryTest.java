package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsCohort;

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

        assertEquals(correctIdT, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, wideSampleId, LimsCohort.WIDE));
        assertEquals(correctIdC, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC, wideSampleId, LimsCohort.WIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId, wideSampleId, LimsCohort.WIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Lims.NOT_AVAILABLE_STRING, wideSampleId, LimsCohort.WIDE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY, wideSampleId, LimsCohort.WIDE));

        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY, coreSampleId, LimsCohort.CORE));
        assertEquals(correctIdT,SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, coreSampleId, LimsCohort.CORE));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId, coreSampleId, LimsCohort.CORE));

        assertNull(correctIdT, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, cpctSampleId, LimsCohort.CPCT));
        assertNull(correctIdC, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC, cpctSampleId, LimsCohort.CPCT));
    }

    @Test
    public void canCheckHospitalPatientId() {
        String coreSampleId = "CORE020000001T";
        String wideSampleId = "WIDE020000001T";

        LimsCohort cohortCORE = LimsCohort.CORE;
        LimsCohort cohortWIDE = LimsCohort.WIDE;

        String hospitalIdNA = Lims.NOT_AVAILABLE_STRING;
        String hospitalIDEmpty = Strings.EMPTY;
        String hospitalId = "1234";

        assertEquals(hospitalIdNA, SampleReportFactory.checkHospitalPatientId(hospitalIdNA, cohortCORE, coreSampleId));
        assertEquals(hospitalIDEmpty, SampleReportFactory.checkHospitalPatientId(hospitalIDEmpty, cohortCORE, coreSampleId));
        assertEquals(hospitalId, SampleReportFactory.checkHospitalPatientId(hospitalId, cohortCORE, coreSampleId));
        assertEquals(hospitalIdNA, SampleReportFactory.checkHospitalPatientId(hospitalIdNA, cohortWIDE, wideSampleId));
    }
}