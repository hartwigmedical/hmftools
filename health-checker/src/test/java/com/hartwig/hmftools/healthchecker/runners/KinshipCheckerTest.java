package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.common.exception.MalformedFileException;
import com.hartwig.hmftools.common.io.dir.RunContext;
import com.hartwig.hmftools.common.io.dir.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.SingleValueResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class KinshipCheckerTest {

    private static final String CORRECT_RUN =
            RunnerTestFunctions.getRunnerResourcePath("kinship") + File.separator + "run";
    private static final String MALFORMED_RUN =
            RunnerTestFunctions.getRunnerResourcePath("kinship") + File.separator + "run2";
    private static final String EMPTY_RUN =
            RunnerTestFunctions.getRunnerResourcePath("kinship") + File.separator + "run3";

    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final double EXPECTED_KINSHIP_VALUE = 0.4748;

    private final KinshipChecker checker = new KinshipChecker();

    @Test
    public void extractDataFromKinship() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(CORRECT_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        Assert.assertEquals(CheckType.KINSHIP, result.getCheckType());
        assertCheck((SingleValueResult) result, Double.toString(EXPECTED_KINSHIP_VALUE));
    }

    @Test
    public void errorYieldsCorrectOutput() {
        final RunContext runContext = TestRunContextFactory.forTest(CORRECT_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        final SingleValueResult result = (SingleValueResult) checker.errorRun(runContext);
        assertNotNull(result.getCheck());
    }

    @Test(expected = MalformedFileException.class)
    public void cannotReadMalformedKinship() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(MALFORMED_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = EmptyFileException.class)
    public void cannotReadFromEmptyKinship() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(EMPTY_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = IOException.class)
    public void cannotReadFromNonExistingKinship() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest("Does not exist", REF_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertCheck(@NotNull final SingleValueResult result, @NotNull final String expectedValue) {
        final HealthCheck healthCheck = result.getCheck();
        assertEquals(expectedValue, healthCheck.getValue());
    }
}
