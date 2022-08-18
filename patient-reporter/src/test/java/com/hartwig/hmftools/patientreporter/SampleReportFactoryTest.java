package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.*;

import org.junit.Test;

public class SampleReportFactoryTest {

    @Test
    public void canTestBiopsyLocation() {
        assertEquals("Ovary", SampleReportFactory.curateBiopsyLocation("Other (please specify below)_ovary"));
        assertEquals("Other", SampleReportFactory.curateBiopsyLocation("Other (please specify below)"));
        assertEquals("Ovary", SampleReportFactory.curateBiopsyLocation("Ovary"));
        assertEquals("Ovary", SampleReportFactory.curateBiopsyLocation("ovary"));
        assertNull(SampleReportFactory.curateBiopsyLocation(null));
    }

    @Test
    public void canTestInterpretReferenceBarcode(){
        assertNull(SampleReportFactory.interpretRefBarcode(null));
        assertEquals("FR123", SampleReportFactory.interpretRefBarcode("FR123"));
        assertEquals("FR123", SampleReportFactory.interpretRefBarcode("FR123"));
        assertEquals("FR123", SampleReportFactory.interpretRefBarcode("FR123-c2f220514"));
        assertEquals("FR123", SampleReportFactory.interpretRefBarcode("FR123_c2f220514"));
    }
}