package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.exception.LineNotFoundException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CoverageCheckerTest {

    private static final String BASE_DIRECTORY = Resources.getResource("coverage").getPath();

    private static final int EXPECTED_NUM_CHECKS = 4;

    private static final String REF_SAMPLE = "sample1";
    private static final double REF_COVERAGE_10X = 0.000037;
    private static final double REF_COVERAGE_20X = 0.00002;
    private static final double REF_COVERAGE_30X = 0.000005;
    private static final double REF_COVERAGE_60X = 0;

    private static final String TUMOR_SAMPLE = "sample2";
    private static final double TUMOR_COVERAGE_10X = 0.000037;
    private static final double TUMOR_COVERAGE_20X = 0.00002;
    private static final double TUMOR_COVERAGE_30X = 0.000005;
    private static final double TUMOR_COVERAGE_60X = 0;

    private static final String EMPTY_SAMPLE = "sample3";
    private static final String INCORRECT_SAMPLE = "sample4";
    private static final String NON_EXISTING_SAMPLE = "sample5";

    private final CoverageChecker checker = new CoverageChecker();

    @Test
    public void correctInputYieldsCorrectOutputForSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.run(runContext);
        assertSomaticResult(result);
    }

    @Test
    public void correctInputYieldsCorrectOutputForSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(BASE_DIRECTORY, REF_SAMPLE);
        final List<HealthCheck> checks = ((MultiValueResult) checker.run(runContext)).getChecks();
        assertEquals(EXPECTED_NUM_CHECKS, checks.size());

        assertCheck(checks, REF_SAMPLE, CoverageCheck.COVERAGE_10X, REF_COVERAGE_10X);
        assertCheck(checks, REF_SAMPLE, CoverageCheck.COVERAGE_20X, REF_COVERAGE_20X);
        assertCheck(checks, REF_SAMPLE, CoverageCheck.COVERAGE_30X, REF_COVERAGE_30X);
        assertCheck(checks, REF_SAMPLE, CoverageCheck.COVERAGE_60X, REF_COVERAGE_60X);
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecksForSomatic() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, result.getTumorSampleChecks().size());
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecksForSingleSample() {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(BASE_DIRECTORY, REF_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
    }

    @Test(expected = EmptyFileException.class)
    public void emptyFileYieldsEmptyFileException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, EMPTY_SAMPLE, EMPTY_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = IOException.class)
    public void nonExistingFileYieldsIOException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, NON_EXISTING_SAMPLE,
                NON_EXISTING_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectRefFileYieldsLineNotFoundException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, INCORRECT_SAMPLE,
                TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectTumorFileYieldsLineNotFoundException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, REF_SAMPLE,
                INCORRECT_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectFilesYieldsLineNotFoundException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(BASE_DIRECTORY, INCORRECT_SAMPLE,
                INCORRECT_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertSomaticResult(@NotNull final BaseResult result) {
        final PatientResult patientResult = (PatientResult) result;

        assertEquals(CheckType.COVERAGE, patientResult.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getTumorSampleChecks().size());

        assertSomaticField(patientResult, CoverageCheck.COVERAGE_10X, REF_COVERAGE_10X, TUMOR_COVERAGE_10X);
        assertSomaticField(patientResult, CoverageCheck.COVERAGE_20X, REF_COVERAGE_20X, TUMOR_COVERAGE_20X);
        assertSomaticField(patientResult, CoverageCheck.COVERAGE_30X, REF_COVERAGE_30X, TUMOR_COVERAGE_30X);
        assertSomaticField(patientResult, CoverageCheck.COVERAGE_60X, REF_COVERAGE_60X, TUMOR_COVERAGE_60X);
    }

    private static void assertSomaticField(@NotNull final PatientResult result, @NotNull final CoverageCheck field,
            final double refValue, final double tumValue) {
        assertCheck(result.getRefSampleChecks(), REF_SAMPLE, field, refValue);
        assertCheck(result.getTumorSampleChecks(), TUMOR_SAMPLE, field, tumValue);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String sampleId,
            @NotNull final CoverageCheck check, final double expectedValue) {
        final Optional<HealthCheck> value = checks.stream().filter(
                p -> p.getCheckName().equals(check.toString())).findFirst();
        assert value.isPresent();

        assertEquals(expectedValue, Double.parseDouble(value.get().getValue()), 1e-10);
        assertEquals(sampleId, value.get().getSampleId());
    }
}
