package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.lims.Lims;

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
        String empty = "";

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
}