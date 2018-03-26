package com.hartwig.hmftools.common.metrics;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.io.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.exception.MalformedFileException;

import org.junit.Test;

public class WGSMetricsFileTest {

    private static final String BASE_DIRECTORY = Resources.getResource("metrics").getPath();
    private static final double EPSILON = 1.0E-10;

    private static final String REF_SAMPLE = "sample1";
    private static final double REF_MEAN_COVERAGE = 0.000856;
    private static final double REF_COVERAGE_10X = 0.000037;
    private static final double REF_COVERAGE_20X = 0.00002;

    private static final String TUMOR_SAMPLE = "sample2";
    private static final double TUMOR_MEAN_COVERAGE = 0.000756;
    private static final double TUMOR_COVERAGE_30X = 0.000005;
    private static final double TUMOR_COVERAGE_60X = 0;

    private static final String EMPTY_SAMPLE = "sample3";
    private static final String INCORRECT_SAMPLE = "sample4";
    private static final String NON_EXISTING_SAMPLE = "sample5";

    @Test
    public void worksForRefAndTumorInput() throws IOException {
        String refFile = WGSMetricsFile.generateFilename(BASE_DIRECTORY, REF_SAMPLE);
        String tumorFile = WGSMetricsFile.generateFilename(BASE_DIRECTORY, TUMOR_SAMPLE);

        WGSMetrics metrics = WGSMetricsFile.read(refFile, tumorFile);

        assertEquals(REF_MEAN_COVERAGE, metrics.refMeanCoverage(), EPSILON);
        assertEquals(REF_COVERAGE_10X, metrics.ref10xCoveragePercentage(), EPSILON);
        assertEquals(REF_COVERAGE_20X, metrics.ref20xCoveragePercentage(), EPSILON);

        Double tumorMeanCov = metrics.tumorMeanCoverage();
        assertNotNull(tumorMeanCov);
        assertEquals(TUMOR_MEAN_COVERAGE, tumorMeanCov, EPSILON);

        Double tumorCov30x = metrics.tumor30xCoveragePercentage();
        assertNotNull(tumorCov30x);
        assertEquals(TUMOR_COVERAGE_30X, tumorCov30x, EPSILON);

        Double tumorCov60x = metrics.tumor60xCoveragePercentage();
        assertNotNull(tumorCov60x);
        assertEquals(TUMOR_COVERAGE_60X, tumorCov60x, EPSILON);
    }

    @Test
    public void worksForRefInputOnly() throws IOException {
        String refFile = WGSMetricsFile.generateFilename(BASE_DIRECTORY, REF_SAMPLE);

        WGSMetrics metrics = WGSMetricsFile.read(refFile);

        assertEquals(REF_MEAN_COVERAGE, metrics.refMeanCoverage(), EPSILON);
        assertEquals(REF_COVERAGE_10X, metrics.ref10xCoveragePercentage(), EPSILON);
        assertEquals(REF_COVERAGE_20X, metrics.ref20xCoveragePercentage(), EPSILON);

        assertNull(metrics.tumorMeanCoverage());
        assertNull(metrics.tumor30xCoveragePercentage());
        assertNull(metrics.tumor60xCoveragePercentage());
    }

    @Test(expected = EmptyFileException.class)
    public void emptyFileYieldsEmptyFileException() throws IOException {
        String file = WGSMetricsFile.generateFilename(BASE_DIRECTORY, EMPTY_SAMPLE);
        WGSMetricsFile.read(file);
    }

    @Test(expected = IOException.class)
    public void nonExistingFileYieldsIOException() throws IOException {
        String file = WGSMetricsFile.generateFilename(BASE_DIRECTORY, NON_EXISTING_SAMPLE);
        WGSMetricsFile.read(file);
    }

    @Test(expected = MalformedFileException.class)
    public void incorrectRefFileYieldsMalformedFileException() throws IOException {
        String file = WGSMetricsFile.generateFilename(BASE_DIRECTORY, INCORRECT_SAMPLE);
        WGSMetricsFile.read(file);
    }
}
