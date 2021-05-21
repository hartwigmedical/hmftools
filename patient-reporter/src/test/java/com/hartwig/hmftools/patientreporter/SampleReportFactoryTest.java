package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.*;

import org.junit.Test;

public class SampleReportFactoryTest {

    @Test
    public void canTestBiopsyLocation() {
        assertEquals("ovary", SampleReportFactory.curateBiopsyLocation("Other (please specify below)_ovary"));
        assertEquals("Other", SampleReportFactory.curateBiopsyLocation("Other (please specify below)"));
        assertEquals("Ovary", SampleReportFactory.curateBiopsyLocation("Ovary"));
        assertNull(SampleReportFactory.curateBiopsyLocation(null));

    }

}