package com.hartwig.hmftools.common.metrics;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class WGSMetricsFileTest {

    private static final String METRICS_DIRECTORY = Resources.getResource("metrics").getPath();
    private static final double EPSILON = 1.0E-10;

    private static final String PV4_FILE = WGSMetricsFile.generateFilename(METRICS_DIRECTORY, "sample.pv4");
    private static final String PV5_FILE = WGSMetricsFile.generateFilename(METRICS_DIRECTORY, "sample.pv5");

    private static final String EMPTY_FILE = WGSMetricsFile.generateFilename(METRICS_DIRECTORY, "sample.empty");
    private static final String INCORRECT_FILE = WGSMetricsFile.generateFilename(METRICS_DIRECTORY, "sample.incorrect");
    private static final String NON_EXISTING_FILE = WGSMetricsFile.generateFilename(METRICS_DIRECTORY, "sample.non-existing");

    @Test
    public void worksForPv4Input() throws IOException {
        WGSMetrics metrics = WGSMetricsFile.read(PV4_FILE);

        assertEquals(0.000856, metrics.meanCoverage(), EPSILON);
        assertEquals(0.257469, metrics.sdCoverage(), EPSILON);
        assertEquals(0, metrics.medianCoverage());
        assertEquals(0, metrics.madCoverage());

        assertNull(metrics.pctExcAdapter());
        assertEquals(0.000585, metrics.pctExcMapQ(), EPSILON);
        assertEquals(0.059484, metrics.pctExcDupe(), EPSILON);
        assertEquals(0.002331, metrics.pctExcUnpaired(), EPSILON);
        assertEquals(0.002378, metrics.pctExcBaseQ(), EPSILON);
        assertEquals(0.020675, metrics.pctExcOverlap(), EPSILON);
        assertEquals(0.001026, metrics.pctExcCapped(), EPSILON);
        assertEquals(0.086479, metrics.pctExcTotal(), EPSILON);

        assertNull(metrics.coverage1xPercentage());
        assertEquals(0.000037, metrics.coverage10xPercentage(), EPSILON);
        assertEquals(0.00002, metrics.coverage20xPercentage(), EPSILON);
        assertEquals(0.000005, metrics.coverage30xPercentage(), EPSILON);
        assertEquals(0D, metrics.coverage60xPercentage(), EPSILON);
    }

    @Test
    public void worksForPv5Input() throws IOException {
        WGSMetrics metrics = WGSMetricsFile.read(PV5_FILE);

        assertEquals(31.440964, metrics.meanCoverage(), EPSILON);
        assertEquals(10.111387, metrics.sdCoverage(), EPSILON);
        assertEquals(32, metrics.medianCoverage());
        assertEquals(5, metrics.madCoverage());

        assertEquals(0.000009, metrics.pctExcAdapter(), EPSILON);
        assertEquals(0.049723, metrics.pctExcMapQ(), EPSILON);
        assertEquals(0.108882, metrics.pctExcDupe(), EPSILON);
        assertEquals(0.000707, metrics.pctExcUnpaired(), EPSILON);
        assertEquals(0.000045, metrics.pctExcBaseQ(), EPSILON);
        assertEquals(0.005271, metrics.pctExcOverlap(), EPSILON);
        assertEquals(0.014476, metrics.pctExcCapped(), EPSILON);
        assertEquals(0.179113, metrics.pctExcTotal(), EPSILON);

        assertEquals(0.99104, metrics.coverage1xPercentage(), EPSILON);
        assertEquals(0.98104, metrics.coverage10xPercentage(), EPSILON);
        assertEquals(0.918671, metrics.coverage20xPercentage(), EPSILON);
        assertEquals(0.630782, metrics.coverage30xPercentage(), EPSILON);
        assertEquals(0.003459, metrics.coverage60xPercentage(), EPSILON);
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
