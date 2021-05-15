package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.junit.Test;

public class PurpleCheckerTest {

    private static final String PURPLE_DIRECTORY = Resources.getResource("purple").getPath();

    @Test
    public void extractDataFromPurpleWorksForSomatic() throws IOException {
        PurpleChecker checker = new PurpleChecker("tumor", PURPLE_DIRECTORY);
        List<QCValue> values = checker.run();

        assertEquals(QCValueType.PURPLE_QC_STATUS, values.get(0).type());
        assertEquals("FAIL_CONTAMINATION", values.get(0).value());
    }
}
