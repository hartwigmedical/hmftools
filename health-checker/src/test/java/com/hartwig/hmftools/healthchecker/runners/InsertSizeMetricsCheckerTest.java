package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.exception.LineNotFoundException;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.context.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.InsertSizeMetricsCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class InsertSizeMetricsCheckerTest {

    private static final String RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("metrics");

    private static final int EXPECTED_NUM_CHECKS = 2;

    private static final String REF_SAMPLE = "sample1";
    private static final String REF_MEDIAN_INSERT_SIZE = "409";
    private static final String REF_WIDTH_OF_70_PERCENT = "195";

    private static final String TUMOR_SAMPLE = "sample2";
    private static final String TUMOR_MEDIAN_INSERT_SIZE = "410";
    private static final String TUMOR_WIDTH_OF_70_PERCENT = "196";

    private static final String EMPTY_SAMPLE = "sample3";
    private static final String INCORRECT_SAMPLE = "sample4";
    private static final String NON_EXISTING_SAMPLE = "sample5";

    private final InsertSizeMetricsChecker checker = new InsertSizeMetricsChecker();

    @Test
    public void correctInputYieldsCorrectOutput() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);
        assertResult(result);
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecks() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, result.getTumorSampleChecks().size());
    }

    @Test(expected = EmptyFileException.class)
    public void emptyFileYieldsEmptyFileException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, EMPTY_SAMPLE, EMPTY_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = IOException.class)
    public void nonExistingFileYieldsIOException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, NON_EXISTING_SAMPLE,
                NON_EXISTING_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectRefFileYieldsLineNotFoundException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, INCORRECT_SAMPLE,
                TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectTumorFileYieldsLineNotFoundException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE,
                INCORRECT_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectFilesYieldsLineNotFoundException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, INCORRECT_SAMPLE,
                INCORRECT_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertResult(@NotNull final BaseResult result) {
        final PatientResult patientResult = (PatientResult) result;
        Assert.assertEquals(CheckType.INSERT_SIZE, patientResult.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getTumorSampleChecks().size());
        assertField(patientResult, InsertSizeMetricsCheck.MAPPING_MEDIAN_INSERT_SIZE.toString(),
                REF_MEDIAN_INSERT_SIZE,
                TUMOR_MEDIAN_INSERT_SIZE);
        assertField(patientResult, InsertSizeMetricsCheck.MAPPING_WIDTH_OF_70_PERCENT.toString(),
                REF_WIDTH_OF_70_PERCENT,
                TUMOR_WIDTH_OF_70_PERCENT);
    }

    private static void assertField(@NotNull final PatientResult result, @NotNull final String field,
            @NotNull final String refValue, @NotNull final String tumValue) {
        assertCheck(result.getRefSampleChecks(), REF_SAMPLE, field, refValue);
        assertCheck(result.getTumorSampleChecks(), TUMOR_SAMPLE, field, tumValue);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String sampleId,
            @NotNull final String checkName, @NotNull final String expectedValue) {
        final Optional<HealthCheck> optCheck = checks.stream().filter(
                p -> p.getCheckName().equals(checkName)).findFirst();
        assert optCheck.isPresent();

        assertEquals(expectedValue, optCheck.get().getValue());
        assertEquals(sampleId, optCheck.get().getSampleId());
    }
}
