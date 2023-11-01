package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.junit.Test;

public class FlagstatCheckerTest {

    private static final String FLAGSTAT_DIRECTORY = Resources.getResource("flagstat").getPath();
    private static final String REF_FLAGSTAT = FLAGSTAT_DIRECTORY + File.separator + "reference.flagstat";
    private static final String TUMOR_FLAGSTAT = FLAGSTAT_DIRECTORY + File.separator + "tumor.flagstat";
    private static final String MALFORMED_FLAGSTAT = FLAGSTAT_DIRECTORY + File.separator + "malformed.flagstat";
    private static final String MISSING_FLAGSTAT = FLAGSTAT_DIRECTORY + File.separator + "doesnotexist.flagstat";

    @Test
    public void extractDataFromFlagstatWorksForSomatic() throws IOException {
        TestFlagstatChecker checker = new TestFlagstatChecker(REF_FLAGSTAT, TUMOR_FLAGSTAT);
        List<QCValue> values = checker.run();

        assertEquals(4, values.size());
        for (QCValue value : values) {
            if (value.Type == QCValueType.REF_PROPORTION_MAPPED) {
                assertEquals("0.993296", value.Value);
            } else if (value.Type == QCValueType.TUM_PROPORTION_MAPPED) {
                assertEquals("0.99683", value.Value);
            } else if (value.Type == QCValueType.REF_PROPORTION_DUPLICATE) {
                assertEquals("0.116625", value.Value);
            } else if (value.Type == QCValueType.TUM_PROPORTION_DUPLICATE) {
                assertEquals("0.163053", value.Value);
            }
        }
    }

    @Test(expected = IOException.class)
    public void malformedYieldsIOException() throws IOException {
        new TestFlagstatChecker(MALFORMED_FLAGSTAT, MALFORMED_FLAGSTAT).run();
    }

    @Test(expected = IOException.class)
    public void missingYieldsIOException() throws IOException {
        new TestFlagstatChecker(MISSING_FLAGSTAT, MISSING_FLAGSTAT).run();
    }
}