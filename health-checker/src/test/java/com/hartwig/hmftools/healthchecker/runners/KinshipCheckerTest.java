package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.exception.MalformedFileException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;
import com.hartwig.hmftools.healthchecker.result.SingleValueResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class KinshipCheckerTest {

    private static final String BASE_DIRECTORY = Resources.getResource("kinship").getPath();

    private static final String CORRECT_RUN = BASE_DIRECTORY + File.separator + "run";
    private static final String MALFORMED_RUN = BASE_DIRECTORY + File.separator + "run2";
    private static final String EMPTY_RUN = BASE_DIRECTORY + File.separator + "run3";

    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final double EXPECTED_KINSHIP_VALUE = 0.4748;

    private final KinshipChecker checker = new KinshipChecker();

    @Test
    public void extractDataFromKinshipWorksForSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(CORRECT_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        Assert.assertEquals(CheckType.KINSHIP, result.getCheckType());
        assertCheck((SingleValueResult) result, Double.toString(EXPECTED_KINSHIP_VALUE));
    }

    @Test
    public void extractDataFromKinshipWorksForSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(CORRECT_RUN, REF_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        assertEquals(CheckType.KINSHIP, result.getCheckType());
        assertTrue(result instanceof NoResult);
    }

    @Test
    public void errorYieldsCorrectOutputForSomatic() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(CORRECT_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        final SingleValueResult result = (SingleValueResult) checker.errorRun(runContext);
        assertNotNull(result.getCheck());
    }

    @Test
    public void errorYieldsCorrectOutputForSingleSample() {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(CORRECT_RUN, REF_SAMPLE);
        final BaseResult result = checker.errorRun(runContext);
        assertTrue(result instanceof NoResult);
    }

    @Test(expected = MalformedFileException.class)
    public void cannotReadMalformedKinship() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(MALFORMED_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = EmptyFileException.class)
    public void cannotReadFromEmptyKinship() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(EMPTY_RUN, REF_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = IOException.class)
    public void cannotReadFromNonExistingKinship() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest("Does not exist", REF_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertCheck(@NotNull final SingleValueResult result, @NotNull final String expectedValue) {
        final HealthCheck healthCheck = result.getCheck();
        assertEquals(expectedValue, healthCheck.getValue());
    }
}
