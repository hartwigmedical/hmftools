package com.hartwig.hmftools.healthcheckreader;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HealthCheckReaderApplicationTest {

    private static final String EXAMPLE_DATA_PATH = Resources.getResource("checks").getPath();
    private static final String EXAMPLE_REPORT_PATH = EXAMPLE_DATA_PATH + File.separator + "example.json";

    @Test
    public void canRunAnalysisOnData() throws IOException {
        new HealthCheckReaderApplication(null, EXAMPLE_DATA_PATH, null).runAnalysis();
    }

    @Test
    public void canRunAnalysisOnReport() throws IOException {
        new HealthCheckReaderApplication(EXAMPLE_REPORT_PATH, null, null).runAnalysis();
    }
}
