package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.exception.LineNotFoundException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.WGSMetricsCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WGSMetricsCheckerTest {

    private static final String RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("metrics");

    private static final int EXPECTED_NUM_CHECKS = 13;

    private static final String REF_SAMPLE = "sample1";
    private static final double REF_COVERAGE_MEAN = 0.000856;
    private static final double REF_COVERAGE_SD = 0.257469;
    private static final double REF_COVERAGE_MEDIAN = 0;
    private static final double REF_COVERAGE_10X = 0.000037;
    private static final double REF_COVERAGE_20X = 0.00002;
    private static final double REF_COVERAGE_30X = 0.000005;
    private static final double REF_COVERAGE_60X = 0;
    private static final double REF_PCT_EXC_MAPQ = 0.000585;
    private static final double REF_PCT_EXC_DUPE = 0.059484;
    private static final double REF_PCT_EXC_UNPAIRED = 0.002331;
    private static final double REF_PCT_EXC_BASEQ = 0.002378;
    private static final double REF_PCT_EXC_OVERLAP = 0.020675;
    private static final double REF_PCT_EXC_TOTAL = 0.086479;

    private static final String TUMOR_SAMPLE = "sample2";
    private static final double TUMOR_COVERAGE_MEAN = 0.000756;
    private static final double TUMOR_COVERAGE_SD = 0.157469;
    private static final double TUMOR_COVERAGE_MEDIAN = 1;
    private static final double TUMOR_COVERAGE_10X = 0.000037;
    private static final double TUMOR_COVERAGE_20X = 0.00002;
    private static final double TUMOR_COVERAGE_30X = 0.000005;
    private static final double TUMOR_COVERAGE_60X = 0;
    private static final double TUMOR_PCT_EXC_MAPQ = 0.000385;
    private static final double TUMOR_PCT_EXC_DUPE = 0.069484;
    private static final double TUMOR_PCT_EXC_UNPAIRED = 0.003331;
    private static final double TUMOR_PCT_EXC_BASEQ = 0.003378;
    private static final double TUMOR_PCT_EXC_OVERLAP = 0.030675;
    private static final double TUMOR_PCT_EXC_TOTAL = 0.096479;

    private static final String EMPTY_SAMPLE = "sample3";
    private static final String INCORRECT_SAMPLE = "sample4";
    private static final String NON_EXISTING_SAMPLE = "sample5";

    private final WGSMetricsChecker checker = new WGSMetricsChecker();

    @Test
    public void correctInputYieldsCorrectOutputForSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.run(runContext);
        assertSomaticResult(result);
    }

    @Test
    public void correctInputYieldsCorrectOutputForSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final List<HealthCheck> checks = ((MultiValueResult) checker.run(runContext)).getChecks();
        assertEquals(EXPECTED_NUM_CHECKS, checks.size());

        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_10X, REF_COVERAGE_10X);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_20X, REF_COVERAGE_20X);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_30X, REF_COVERAGE_30X);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_60X, REF_COVERAGE_60X);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_MEAN, REF_COVERAGE_MEAN);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_MEDIAN, REF_COVERAGE_MEDIAN);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_PCT_EXC_BASEQ, REF_PCT_EXC_BASEQ);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_PCT_EXC_DUPE, REF_PCT_EXC_DUPE);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_PCT_EXC_OVERLAP, REF_PCT_EXC_OVERLAP);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_PCT_EXC_MAPQ, REF_PCT_EXC_MAPQ);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_PCT_EXC_TOTAL, REF_PCT_EXC_TOTAL);
        assertCheck(checks, REF_SAMPLE, WGSMetricsCheck.COVERAGE_PCT_EXC_UNPAIRED, REF_PCT_EXC_UNPAIRED);
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecksForSomatic() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, result.getTumorSampleChecks().size());
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecksForSingleSample() {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
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

    private static void assertSomaticResult(@NotNull final BaseResult result) {
        final PatientResult patientResult = (PatientResult) result;

        assertEquals(CheckType.WGS_METRICS, patientResult.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getTumorSampleChecks().size());
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_MEAN, REF_COVERAGE_MEAN, TUMOR_COVERAGE_MEAN);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_MEDIAN, REF_COVERAGE_MEDIAN, TUMOR_COVERAGE_MEDIAN);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_SD, REF_COVERAGE_SD, TUMOR_COVERAGE_SD);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_10X, REF_COVERAGE_10X, TUMOR_COVERAGE_10X);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_20X, REF_COVERAGE_20X, TUMOR_COVERAGE_20X);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_30X, REF_COVERAGE_30X, TUMOR_COVERAGE_30X);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_60X, REF_COVERAGE_60X, TUMOR_COVERAGE_60X);

        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_BASEQ, REF_PCT_EXC_BASEQ,
                TUMOR_PCT_EXC_BASEQ);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_DUPE, REF_PCT_EXC_DUPE, TUMOR_PCT_EXC_DUPE);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_MAPQ, REF_PCT_EXC_MAPQ, TUMOR_PCT_EXC_MAPQ);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_OVERLAP, REF_PCT_EXC_OVERLAP,
                TUMOR_PCT_EXC_OVERLAP);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_UNPAIRED, REF_PCT_EXC_UNPAIRED,
                TUMOR_PCT_EXC_UNPAIRED);
        assertSomaticField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_TOTAL, REF_PCT_EXC_TOTAL,
                TUMOR_PCT_EXC_TOTAL);
    }

    private static void assertSomaticField(@NotNull final PatientResult result, @NotNull final WGSMetricsCheck field,
            final double refValue, final double tumValue) {
        assertCheck(result.getRefSampleChecks(), REF_SAMPLE, field, refValue);
        assertCheck(result.getTumorSampleChecks(), TUMOR_SAMPLE, field, tumValue);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String sampleId,
            @NotNull final WGSMetricsCheck check, final double expectedValue) {
        final Optional<HealthCheck> value = checks.stream().filter(
                p -> p.getCheckName().equals(check.toString())).findFirst();
        assert value.isPresent();

        assertEquals(expectedValue, Double.parseDouble(value.get().getValue()), 1e-10);
        assertEquals(sampleId, value.get().getSampleId());
    }
}
