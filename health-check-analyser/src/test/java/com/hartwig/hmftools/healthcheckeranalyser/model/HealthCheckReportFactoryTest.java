package com.hartwig.hmftools.healthcheckeranalyser.model;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HealthCheckReportFactoryTest {

    private static final String EXAMPLE_REPORT_PATH =
            Resources.getResource("checks").getPath() + File.separator + "example.json";

    @Test
    public void canConvertReportToJavaModel() throws IOException {
        HealthCheckReport report = HealthCheckReportFactory.fromHealthCheckReport(EXAMPLE_REPORT_PATH);
        assertNotNull(report);
    }
}