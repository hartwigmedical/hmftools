package com.hartwig.hmftools.healthcheckeranalyser;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HealthCheckerAnalysisApplicationTest {

    private static final String EXAMPLE_REPORT_PATH =
            Resources.getResource("checks").getPath() + File.separator + "example.json";

    @Test
    public void canRunAnalysis() throws IOException {
        new HealthCheckerAnalysisApplication(EXAMPLE_REPORT_PATH).runAnalysis();
    }
}
