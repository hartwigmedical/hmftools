package com.hartwig.hmftools.healthcheckeranalyser;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HealthCheckerAnalysisApplicationTest {

    private static final String EXAMPLE_REPORT_PATH = Resources.getResource("checks").getPath();

    @Test
    public void canRunAnalysis() throws IOException {
        String csvOut = "/Users/kduyvesteyn/hmf/tmp/checks.csv";
        new HealthCheckerAnalysisApplication(EXAMPLE_REPORT_PATH, csvOut).runAnalysis();
    }
}
