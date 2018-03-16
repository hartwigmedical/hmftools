package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.MalformedFileException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.SingleValueResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class AmberCheckerTest {
    private static final String BASE_DIRECTORY = Resources.getResource("").getPath();
    private static final String REF_SAMPLE = "refSample";
    private static final String TUMOR_SAMPLE = "tumorSample";

    private final AmberChecker checker = new AmberChecker();

    @Test
    public void extractDataFromAmberWorksForSomatic() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.run(runContext);

        Assert.assertEquals(CheckType.AMBER, result.getCheckType());
        final HealthCheck check = ((SingleValueResult) result).getCheck();
        assertCheck(check, "0.4951");
    }

    @Test(expected = MalformedFileException.class)
    public void testMalformed() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, "malformed");
        checker.run(runContext);
    }

    @Test(expected = IOException.class)
    public void testMissing() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, "missing");
        checker.run(runContext);
    }

    private static void assertCheck(@NotNull final HealthCheck check, final String expectedValue) {
        assertEquals(expectedValue, check.getValue());
        assertEquals(AmberCheck.MEAN_BAF.toString(), check.getCheckName());
    }
}
