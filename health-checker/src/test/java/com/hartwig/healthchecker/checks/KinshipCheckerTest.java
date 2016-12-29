package com.hartwig.healthchecker.checks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.healthchecker.common.checks.CheckType;
import com.hartwig.healthchecker.common.checks.HealthCheck;
import com.hartwig.healthchecker.common.exception.EmptyFileException;
import com.hartwig.healthchecker.common.exception.HealthChecksException;
import com.hartwig.healthchecker.common.exception.MalformedFileException;
import com.hartwig.healthchecker.common.io.dir.RunContext;
import com.hartwig.healthchecker.common.io.dir.TestRunContextFactory;
import com.hartwig.healthchecker.common.result.BaseResult;
import com.hartwig.healthchecker.common.result.SingleValueResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class KinshipCheckerTest {

    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final String CORRECT_RUN = Resources.getResource("checks/kinship/run").getPath();
    private static final String MALFORMED_RUN = Resources.getResource("checks/kinship/run2").getPath();
    private static final String EMPTY_RUN = Resources.getResource("checks/kinship/run3").getPath();

    private static final double EXPECTED_KINSHIP_VALUE = 0.4748;

    private final KinshipChecker checker = new KinshipChecker();

    @Test
    public void extractDataFromKinship() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(CORRECT_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        assertEquals(CheckType.KINSHIP, result.getCheckType());
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

    private static void assertCheck(@NotNull final SingleValueResult result,
            @NotNull final String expectedValue) {
        final HealthCheck healthCheck = result.getCheck();
        assertEquals(expectedValue, healthCheck.getValue());
    }
}
