package com.hartwig.hmftools.healthchecker;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.runners.HealthCheckSampleConfiguration;

import org.junit.Test;

public class HealthChecksApplicationTest {

    private static final String BASE_DIR = Resources.getResource("").getPath();

    private static final String OUTPUT_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final boolean WRITE_EVALUATION_FILE = false;

    @Test
    public void runHealthCheckerInSomaticMode() throws IOException {
        HealthChecksApplication app = new HealthChecksApplication(HealthCheckSampleConfiguration.of("reference",
                BASE_DIR + "metrics" + "/reference.wgsmetrics",
                BASE_DIR + "flagstat/reference.flagstat"),
                HealthCheckSampleConfiguration.of("tumor", BASE_DIR + "metrics/tumor.wgsmetrics", BASE_DIR + "flagstat/tumor.flagstat"),
                BASE_DIR + "purple",
                OUTPUT_DIR);

        app.run(WRITE_EVALUATION_FILE);
    }

    @Test
    public void runHealthCheckerInGermlineOnlyMode() throws IOException {
        HealthChecksApplication app = new HealthChecksApplication(HealthCheckSampleConfiguration.of("reference",
                BASE_DIR + "metrics" + "/reference.wgsmetrics",
                BASE_DIR + "flagstat/reference.flagstat"), null, null, OUTPUT_DIR);

        app.run(WRITE_EVALUATION_FILE);
    }

    @Test
    public void runHealthCheckerInTumorOnlyMode() throws IOException {
        HealthChecksApplication app = new HealthChecksApplication(null,
                HealthCheckSampleConfiguration.of("tumor", BASE_DIR + "metrics/tumor.wgsmetrics", BASE_DIR + "flagstat/tumor.flagstat"),
                BASE_DIR + "purple",
                OUTPUT_DIR);

        app.run(WRITE_EVALUATION_FILE);
    }
}