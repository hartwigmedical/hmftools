package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.junit.Test;

public class FlagstatCheckerTest {

    private static final String FLAGSTAT_DIRECTORY = Resources.getResource("flagstat").getPath();
    private static final String REF_FLAGSTAT = FLAGSTAT_DIRECTORY + "/reference.flagstat";
    private static final String TUM_FLAGSTAT = FLAGSTAT_DIRECTORY + "/tumor.flagstat";
    private static final String MALFORMED_FLAGSTAT = FLAGSTAT_DIRECTORY + "/malformed.flagstat";

    @Test
    public void extractDataFromFlagstatWorksForSomatic() throws IOException {
        FlagstatChecker checker = new FlagstatChecker(REF_FLAGSTAT, TUM_FLAGSTAT);
        List<QCValue> values = checker.run();

        assertEquals(2, values.size());
        for (QCValue value : values) {
            if (value.type() == QCValueType.REF_PROPORTION_MAPPED) {
                assertEquals("0.9932958382082092", value.value());
            } else if (value.type() == QCValueType.TUM_PROPORTION_MAPPED) {
                assertEquals("0.9968295130748326", value.value());
            }
        }
    }

    @Test(expected = IOException.class)
    public void malformedYieldsIOException() throws IOException {
        FlagstatChecker checker = new FlagstatChecker(MALFORMED_FLAGSTAT, MALFORMED_FLAGSTAT);
        checker.run();
    }

    @Test(expected = IOException.class)
    public void missingYieldsIOException() throws IOException {
        FlagstatChecker checker = new FlagstatChecker("missing", "missing");
        checker.run();
    }
}