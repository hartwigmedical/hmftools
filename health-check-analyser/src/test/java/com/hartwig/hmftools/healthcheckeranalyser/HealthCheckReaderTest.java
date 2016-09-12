package com.hartwig.hmftools.healthcheckeranalyser;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReport;

import org.junit.Test;

public class HealthCheckReaderTest {

    private static final String RUN_DIRECTORY = Resources.getResource("checks").getPath();
    private static final String EXAMPLE_REPORT = RUN_DIRECTORY + File.separator + "example.json";

    @Test
    public void canConvertReportToJavaModel() throws IOException {
        HealthCheckReport report = HealthCheckReader.readHealthCheckOutput(EXAMPLE_REPORT);
        assertNotNull(report);
    }
}