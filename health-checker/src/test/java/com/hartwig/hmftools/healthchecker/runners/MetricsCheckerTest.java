package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.junit.Test;

public class MetricsCheckerTest {

    private static final String METRICS_DIRECTORY = Resources.getResource("metrics").getPath();

    @Test
    public void extractDataFromMetricsWorksForSomatic() throws IOException {
        MetricsChecker checker = new MetricsChecker("reference", "tumor", METRICS_DIRECTORY);
        List<QCValue> values = checker.run();

        assertEquals(4, values.size());
        for (QCValue value : values) {
            if (value.type() == QCValueType.REF_COVERAGE_10X) {
                assertEquals("0.98261", value.value());
            } else if (value.type() == QCValueType.REF_COVERAGE_20X) {
                assertEquals("0.980701", value.value());
            } else if (value.type() == QCValueType.TUMOR_COVERAGE_30X) {
                assertEquals("0.978779", value.value());
            } else if (value.type() == QCValueType.TUMOR_COVERAGE_60X) {
                assertEquals("0.950264", value.value());
            }
        }
    }

    @Test
    public void extractDataFromMetricsWorksForSingleSample() throws IOException {
        MetricsChecker checker = new MetricsChecker("reference", null, METRICS_DIRECTORY);
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
        MetricsChecker checker = new MetricsChecker("malformed", null, METRICS_DIRECTORY);
        checker.run();
    }

    @Test(expected = IOException.class)
    public void missingYieldsIOException() throws IOException {
        MetricsChecker checker = new MetricsChecker("missing", null, METRICS_DIRECTORY);
        checker.run();
    }
}