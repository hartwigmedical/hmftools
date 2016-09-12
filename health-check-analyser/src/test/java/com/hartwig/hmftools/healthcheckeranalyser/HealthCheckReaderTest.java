package com.hartwig.hmftools.healthcheckeranalyser;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReport;

import org.junit.Test;

public class HealthCheckReaderTest {

    private static final String EXAMPLE_REPORT_PATH =
            Resources.getResource("checks").getPath() + File.separator + "example.json";

    @Test
    public void canConvertReportToJavaModel() throws IOException {
        HealthCheckReport report = HealthCheckReader.readHealthCheckOutput(EXAMPLE_REPORT_PATH);
        assertNotNull(report);
    }
}