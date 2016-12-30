package com.hartwig.hmftools.healthchecker.report;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.hartwig.hmftools.healthchecker.exception.HealthChecksException;

import org.junit.Test;

public class ReportsFlyweightTest {

    private static final String BLA = "BLA";
    private static final String STDOUT = "stdout";
    private static final String JSON = "json";

    @Test
    public void getReport() throws IOException, HealthChecksException {
        Report report = ReportsFlyweight.getInstance().getReport(STDOUT);
        assertNotNull(report);
        assertTrue(report instanceof StandardOutputReport);

        report = ReportsFlyweight.getInstance().getReport(JSON);
        assertNotNull(report);
        assertTrue(report instanceof JsonReport);

        report = ReportsFlyweight.getInstance().getReport(BLA);
        assertNotNull(report);
        assertTrue(report instanceof StandardOutputReport);
    }
}
