package com.hartwig.hmftools.healthcheckreader;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HealthCheckReaderApplicationTest {

    private static final String EXAMPLE_REPORT_PATH =
            Resources.getResource("checks").getPath() + File.separator + "example.json";
    private static final String EXAMPLE_DATA_PATH = Resources.getResource("checks").getPath();

    @Test
    public void canRunAnalysisOnData() throws IOException {
        String csvOut = "/Users/kduyvesteyn/hmf/tmp/checks.csv";
        new HealthCheckReaderApplication(null, EXAMPLE_DATA_PATH, csvOut).runAnalysis();
    }

    @Test
    public void canRunAnalysisOnReport() throws IOException {
        new HealthCheckReaderApplication(EXAMPLE_REPORT_PATH, null, null).runAnalysis();
    }
}
