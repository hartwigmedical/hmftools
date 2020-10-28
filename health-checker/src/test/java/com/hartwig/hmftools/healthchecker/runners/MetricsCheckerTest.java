package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.junit.Test;

public class MetricsCheckerTest {

    private static final String METRICS_DIRECTORY = Resources.getResource("metrics").getPath();

    private static final String REF_WGS_METRICS_FILE = METRICS_DIRECTORY + File.separator + "reference.wgsmetrics";
    private static final String TUM_WGS_METRICS_FILE = METRICS_DIRECTORY + File.separator + "tumor.wgsmetrics";
    private static final String MALFORMED_WGS_METRICS_FILE = METRICS_DIRECTORY + File.separator + "malformed.wgsmetrics";
    private static final String MISSING_WGS_METRICS_FILE = METRICS_DIRECTORY + File.separator + "doesnotexist.wgsmetrics";

    @Test
    public void extractDataFromMetricsWorksForSomatic() throws IOException {
        MetricsChecker checker = new MetricsChecker(REF_WGS_METRICS_FILE, TUM_WGS_METRICS_FILE);
        List<QCValue> values = checker.run();

        assertEquals(4, values.size());
        for (QCValue value : values) {
            if (value.type() == QCValueType.REF_COVERAGE_10X) {
                assertEquals("0.98261", value.value());
            } else if (value.type() == QCValueType.REF_COVERAGE_20X) {
                assertEquals("0.980701", value.value());
            } else if (value.type() == QCValueType.TUM_COVERAGE_30X) {
                assertEquals("0.978779", value.value());
            } else if (value.type() == QCValueType.TUM_COVERAGE_60X) {
                assertEquals("0.950264", value.value());
            }
        }
    }

    @Test
    public void extractDataFromMetricsWorksForSingleSample() throws IOException {
        MetricsChecker checker = new MetricsChecker(REF_WGS_METRICS_FILE, null);
        List<QCValue> values = checker.run();

        assertEquals(2, values.size());
        for (QCValue value : values) {
            if (value.type() == QCValueType.REF_COVERAGE_10X) {
                assertEquals("0.98261", value.value());
            } else if (value.type() == QCValueType.REF_COVERAGE_20X) {
                assertEquals("0.980701", value.value());
            }
        }
    }

    @Test(expected = IOException.class)
    public void malformedYieldsIOException() throws IOException {
        new MetricsChecker(MALFORMED_WGS_METRICS_FILE, MALFORMED_WGS_METRICS_FILE).run();
    }

    @Test(expected = IOException.class)
    public void missingYieldsIOException() throws IOException {
        new MetricsChecker(MISSING_WGS_METRICS_FILE, MISSING_WGS_METRICS_FILE).run();
    }
}