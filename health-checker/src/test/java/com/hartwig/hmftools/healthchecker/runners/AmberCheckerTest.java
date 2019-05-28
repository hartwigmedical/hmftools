package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.junit.Test;

public class AmberCheckerTest {
    private static final String AMBER_DIRECTORY = Resources.getResource("amber").getPath();

    @Test
    public void extractDataFromAmberWorksForSomatic() throws IOException {
        final AmberChecker checker = new AmberChecker("tumor", AMBER_DIRECTORY);
        final List<QCValue> values = checker.run();

        for (QCValue value : values) {
            if (value.type() == QCValueType.AMBER_MEAN_BAF) {
                assertEquals("0.4951", value.value());
            } else if (value.type() == QCValueType.AMBER_CONTAMINATION) {
                assertEquals("0.001", value.value());
            }
        }
    }

    @Test(expected = IOException.class)
    public void malformedYieldsIOException() throws IOException {
        final AmberChecker checker = new AmberChecker("malformed", AMBER_DIRECTORY);
        checker.run();
    }

    @Test(expected = IOException.class)
    public void missingYieldsIOException() throws IOException {
        final AmberChecker checker = new AmberChecker("missing", AMBER_DIRECTORY);
        checker.run();
    }
}
