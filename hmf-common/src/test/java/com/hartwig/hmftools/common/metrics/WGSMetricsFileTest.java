package com.hartwig.hmftools.common.metrics;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class WGSMetricsFileTest {

    private static final String METRICS_DIRECTORY = Resources.getResource("metrics").getPath();
    private static final double EPSILON = 1.0E-10;

    private static final String PV4_REF_SAMPLE = "sample1_pv4";
    private static final String PV4_REF_FILE = METRICS_DIRECTORY + File.separator + PV4_REF_SAMPLE + ".wgsmetrics";
    private static final double PV4_REF_MEAN_COVERAGE = 0.000856;
    private static final double PV4_REF_COVERAGE_10X = 0.000037;
    private static final double PV4_REF_COVERAGE_20X = 0.00002;

    private static final String PV4_TUMOR_SAMPLE = "sample2_pv4";
    private static final String PV4_TUMOR_FILE = METRICS_DIRECTORY + File.separator + PV4_TUMOR_SAMPLE + ".wgsmetrics";
    private static final double PV4_TUMOR_MEAN_COVERAGE = 0.000756;
    private static final double PV4_TUMOR_COVERAGE_30X = 0.000005;
    private static final double PV4_TUMOR_COVERAGE_60X = 0;

    private static final String PV5_REF_SAMPLE = "sample1_pv5";
    private static final String PV5_REF_FILE = METRICS_DIRECTORY + File.separator + PV5_REF_SAMPLE + ".wgsmetrics";
    private static final double PV5_REF_MEAN_COVERAGE = 37.190252;
    private static final double PV5_REF_COVERAGE_10X = 0.979308;
    private static final double PV5_REF_COVERAGE_20X = 0.96964;

    private static final String PV5_TUMOR_SAMPLE = "sample2_pv5";
    private static final String PV5_TUMOR_FILE = METRICS_DIRECTORY + File.separator + PV5_TUMOR_SAMPLE + ".wgsmetrics";
    private static final double PV5_TUMOR_MEAN_COVERAGE = 122.955213;
    private static final double PV5_TUMOR_COVERAGE_30X = 0.980146;
    private static final double PV5_TUMOR_COVERAGE_60X = 0.974498;

    private static final String EMPTY_FILE = METRICS_DIRECTORY + File.separator + "sample3" + ".wgsmetrics";
    private static final String INCORRECT_FILE = METRICS_DIRECTORY + File.separator + "sample4" + ".wgsmetrics";
    private static final String NON_EXISTING_FILE = METRICS_DIRECTORY + File.separator + "sample5" + ".wgsmetrics";

    @Test
    public void worksForPv4RefAndTumorInput() throws IOException {
        WGSMetrics metrics = WGSMetricsFile.read(PV4_REF_FILE, PV4_TUMOR_FILE);

        assertEquals(PV4_REF_MEAN_COVERAGE, metrics.refMeanCoverage(), EPSILON);
        assertEquals(PV4_REF_COVERAGE_10X, metrics.ref10xCoveragePercentage(), EPSILON);
        assertEquals(PV4_REF_COVERAGE_20X, metrics.ref20xCoveragePercentage(), EPSILON);

        Double tumorMeanCov = metrics.tumorMeanCoverage();
        assertNotNull(tumorMeanCov);
        assertEquals(PV4_TUMOR_MEAN_COVERAGE, tumorMeanCov, EPSILON);

        Double tumorCov30x = metrics.tumor30xCoveragePercentage();
        assertNotNull(tumorCov30x);
        assertEquals(PV4_TUMOR_COVERAGE_30X, tumorCov30x, EPSILON);

        Double tumorCov60x = metrics.tumor60xCoveragePercentage();
        assertNotNull(tumorCov60x);
        assertEquals(PV4_TUMOR_COVERAGE_60X, tumorCov60x, EPSILON);
    }

    @Test
    public void worksForPv4RefInputOnly() throws IOException {
        WGSMetrics metrics = WGSMetricsFile.read(PV4_REF_FILE);

        assertEquals(PV4_REF_MEAN_COVERAGE, metrics.refMeanCoverage(), EPSILON);
        assertEquals(PV4_REF_COVERAGE_10X, metrics.ref10xCoveragePercentage(), EPSILON);
        assertEquals(PV4_REF_COVERAGE_20X, metrics.ref20xCoveragePercentage(), EPSILON);

        assertNull(metrics.tumorMeanCoverage());
        assertNull(metrics.tumor30xCoveragePercentage());
        assertNull(metrics.tumor60xCoveragePercentage());
    }

    @Test
    public void worksForPv5RefAndTumorInput() throws IOException {
        WGSMetrics metrics = WGSMetricsFile.read(PV5_REF_FILE,PV5_TUMOR_FILE);

        assertEquals(PV5_REF_MEAN_COVERAGE, metrics.refMeanCoverage(), EPSILON);
        assertEquals(PV5_REF_COVERAGE_10X, metrics.ref10xCoveragePercentage(), EPSILON);
        assertEquals(PV5_REF_COVERAGE_20X, metrics.ref20xCoveragePercentage(), EPSILON);

        Double tumorMeanCov = metrics.tumorMeanCoverage();
        assertNotNull(tumorMeanCov);
        assertEquals(PV5_TUMOR_MEAN_COVERAGE, tumorMeanCov, EPSILON);

        Double tumorCov30x = metrics.tumor30xCoveragePercentage();
        assertNotNull(tumorCov30x);
        assertEquals(PV5_TUMOR_COVERAGE_30X, tumorCov30x, EPSILON);

        Double tumorCov60x = metrics.tumor60xCoveragePercentage();
        assertNotNull(tumorCov60x);
        assertEquals(PV5_TUMOR_COVERAGE_60X, tumorCov60x, EPSILON);
    }

    @Test(expected = IOException.class)
    public void emptyFileYieldsIOException() throws IOException {
        WGSMetricsFile.read(EMPTY_FILE);
    }

    @Test(expected = IOException.class)
    public void nonExistingFileYieldsIOException() throws IOException {
        WGSMetricsFile.read(NON_EXISTING_FILE);
    }

    @Test(expected = IOException.class)
    public void incorrectRefFileYieldsIOException() throws IOException {
        WGSMetricsFile.read(INCORRECT_FILE);
    }
}
