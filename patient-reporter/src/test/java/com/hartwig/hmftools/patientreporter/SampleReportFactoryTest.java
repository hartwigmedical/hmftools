package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.*;

import org.junit.Test;

public class SampleReportFactoryTest {

    @Test
    public void canConvertPaTissueId() {
        String sampleID = "WIDE020000001T";
        String correctIdT = "T20-72346";
        String correctIdC = "C18-00124";
        String wrongId = "BGr-12111";
        String NAId = "N/A";
        String emptyTissueId = "";

        assertTrue(SampleReportFactory.reportHospitalTissueIdPA(correctIdT, sampleID));
        assertTrue(SampleReportFactory.reportHospitalTissueIdPA(correctIdC, sampleID));
        assertFalse(SampleReportFactory.reportHospitalTissueIdPA(wrongId, sampleID));
        assertFalse(SampleReportFactory.reportHospitalTissueIdPA(NAId, sampleID));
        assertFalse(SampleReportFactory.reportHospitalTissueIdPA(emptyTissueId, sampleID));

    }

}