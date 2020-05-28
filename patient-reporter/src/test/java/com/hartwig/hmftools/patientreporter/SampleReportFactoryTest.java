package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsStudy;

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

        assertEquals(correctIdT, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, wideSampleId));
        assertEquals(correctIdC, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC, wideSampleId));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId, wideSampleId));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Lims.NOT_AVAILABLE_STRING, wideSampleId));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY, wideSampleId));

        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(Strings.EMPTY, coreSampleId));
        assertEquals(correctIdT,SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, coreSampleId));
        assertNull(SampleReportFactory.toHospitalPathologySampleIdForReport(wrongId, coreSampleId));

        assertNull(correctIdT, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdT, cpctSampleId));
        assertNull(correctIdC, SampleReportFactory.toHospitalPathologySampleIdForReport(correctIdC, cpctSampleId));
    }

    @Test
    public void canCheckHospitalPatientId() {
        String coreSampleId = "CORE020000001T";
        String wideSampleId = "WIDE020000001T";

        LimsStudy typeCORE = LimsStudy.CORE;
        LimsStudy typeWIDE = LimsStudy.WIDE;

        String hositalIdNA = Lims.NOT_AVAILABLE_STRING;
        String hospitalIDempty = Strings.EMPTY;
        String hospitalId = "1234";

        assertEquals(hositalIdNA, SampleReportFactory.checkHospitalPatientId(hositalIdNA, typeCORE, coreSampleId));
        assertEquals(hospitalIDempty, SampleReportFactory.checkHospitalPatientId(hospitalIDempty, typeCORE, coreSampleId));
        assertEquals(hospitalId, SampleReportFactory.checkHospitalPatientId(hospitalId, typeCORE, coreSampleId));
        assertEquals(hositalIdNA, SampleReportFactory.checkHospitalPatientId(hositalIdNA, typeWIDE, wideSampleId));

    }
}