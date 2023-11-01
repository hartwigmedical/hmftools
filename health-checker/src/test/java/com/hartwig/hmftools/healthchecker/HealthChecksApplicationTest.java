package com.hartwig.hmftools.healthchecker;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.runners.HealthCheckSampleConfiguration;

import org.junit.Test;

public class HealthChecksApplicationTest
{
    private static final String BASE_DIR = Resources.getResource("").getPath();

    private static final String OUTPUT_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    private static final HealthCheckSampleConfiguration REF_CONFIG = new HealthCheckSampleConfiguration(
            "reference",
            BASE_DIR + "metrics/reference.wgsmetrics",
            BASE_DIR + "flagstat/reference.flagstat");

    private static final HealthCheckSampleConfiguration TUMOR_CONFIG = new HealthCheckSampleConfiguration(
            "tumor",
            BASE_DIR + "metrics/tumor.wgsmetrics",
            BASE_DIR + "flagstat/tumor.flagstat");

    @Test
    public void runHealthCheckerInSomaticMode() throws IOException
    {
        HealthChecksApplication app = new HealthChecksApplication(REF_CONFIG, TUMOR_CONFIG, BASE_DIR + "purple", OUTPUT_DIR);
        app.run();
    }

    @Test
    public void runHealthCheckerInGermlineOnlyMode() throws IOException
    {
        HealthChecksApplication app = new HealthChecksApplication(REF_CONFIG, null, null, OUTPUT_DIR);
        app.run();
    }

    @Test
    public void runHealthCheckerInTumorOnlyMode() throws IOException
    {
        HealthChecksApplication app = new HealthChecksApplication(null, TUMOR_CONFIG, BASE_DIR + "purple", OUTPUT_DIR);
        app.run();
    }
}