package com.hartwig.hmftools.common.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.healthchecker.result.QCValue;
import com.hartwig.hmftools.common.healthchecker.result.QCValueType;

import org.junit.Test;

public class PurpleCheckerTest {

    private static final String PURPLE_DIRECTORY = Resources.getResource("healthchecker/purple").getPath();

    @Test
    public void extractDataFromPurpleWorksForSomatic() throws IOException {
        final PurpleChecker checker = new PurpleChecker("tumor", PURPLE_DIRECTORY);
        final List<QCValue> values = checker.run();

        assertEquals(QCValueType.PURPLE_QC_STATUS, values.get(0).type());
        assertEquals("FAIL_GENDER", values.get(0).value());
    }

    @Test (expected = IOException.class)
    public void malformedYieldsIOException() throws IOException {
        final PurpleChecker checker = new PurpleChecker("malformed", PURPLE_DIRECTORY);
        checker.run();
    }

    @Test (expected = IOException.class)
    public void missingYieldsIOException() throws IOException {
        final PurpleChecker checker = new PurpleChecker("missing", PURPLE_DIRECTORY);
        checker.run();
    }
}
