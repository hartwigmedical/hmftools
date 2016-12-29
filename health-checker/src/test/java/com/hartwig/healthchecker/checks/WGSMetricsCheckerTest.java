package com.hartwig.healthchecker.checks;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.healthchecker.common.checks.CheckType;
import com.hartwig.healthchecker.common.checks.HealthCheck;
import com.hartwig.healthchecker.common.exception.EmptyFileException;
import com.hartwig.healthchecker.common.exception.HealthChecksException;
import com.hartwig.healthchecker.common.exception.LineNotFoundException;
import com.hartwig.healthchecker.common.io.dir.RunContext;
import com.hartwig.healthchecker.common.io.dir.TestRunContextFactory;
import com.hartwig.healthchecker.common.result.BaseResult;
import com.hartwig.healthchecker.common.result.PatientResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WGSMetricsCheckerTest {

    private static final String RUN_DIRECTORY = Resources.getResource("checks/metrics").getPath();

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
    public void correctInputYieldsCorrectOutput() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.run(runContext);
        assertResult(result);
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecks() {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, result.getTumorSampleChecks().size());
    }

    @Test(expected = EmptyFileException.class)
    public void emptyFileYieldsEmptyFileException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, EMPTY_SAMPLE, EMPTY_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = IOException.class)
    public void nonExistingFileYieldsIOException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, NON_EXISTING_SAMPLE,
                NON_EXISTING_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectRefFileYieldsLineNotFoundException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, INCORRECT_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectTumorFileYieldsLineNotFoundException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, INCORRECT_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = LineNotFoundException.class)
    public void incorrectFilesYieldsLineNotFoundException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, INCORRECT_SAMPLE, INCORRECT_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertResult(@NotNull final BaseResult result) {
        final PatientResult patientResult = (PatientResult) result;

        assertEquals(CheckType.WGS_METRICS, patientResult.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getTumorSampleChecks().size());
        assertField(patientResult, WGSMetricsCheck.COVERAGE_MEAN.name(), REF_COVERAGE_MEAN, TUMOR_COVERAGE_MEAN);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_MEDIAN.name(), REF_COVERAGE_MEDIAN, TUMOR_COVERAGE_MEDIAN);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_SD.name(), REF_COVERAGE_SD, TUMOR_COVERAGE_SD);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_10X.name(), REF_COVERAGE_10X, TUMOR_COVERAGE_10X);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_20X.name(), REF_COVERAGE_20X, TUMOR_COVERAGE_20X);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_30X.name(), REF_COVERAGE_30X, TUMOR_COVERAGE_30X);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_60X.name(), REF_COVERAGE_60X, TUMOR_COVERAGE_60X);

        assertField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_BASEQ.name(), REF_PCT_EXC_BASEQ,
                TUMOR_PCT_EXC_BASEQ);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_DUPE.name(), REF_PCT_EXC_DUPE, TUMOR_PCT_EXC_DUPE);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_MAPQ.name(), REF_PCT_EXC_MAPQ, TUMOR_PCT_EXC_MAPQ);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_OVERLAP.name(), REF_PCT_EXC_OVERLAP,
                TUMOR_PCT_EXC_OVERLAP);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_UNPAIRED.name(), REF_PCT_EXC_UNPAIRED,
                TUMOR_PCT_EXC_UNPAIRED);
        assertField(patientResult, WGSMetricsCheck.COVERAGE_PCT_EXC_TOTAL.name(), REF_PCT_EXC_TOTAL,
                TUMOR_PCT_EXC_TOTAL);
    }

    private static void assertField(@NotNull final PatientResult result, @NotNull final String field,
            final double refValue, final double tumValue) {
        assertCheck(result.getRefSampleChecks(), REF_SAMPLE, field, refValue);
        assertCheck(result.getTumorSampleChecks(), TUMOR_SAMPLE, field, tumValue);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String sampleId,
            @NotNull final String checkName, final double expectedValue) {
        final Optional<HealthCheck> value = checks.stream().filter(
                p -> p.getCheckName().equals(checkName)).findFirst();
        assert value.isPresent();

        assertEquals(expectedValue, Double.parseDouble(value.get().getValue()), 1e-10);
        assertEquals(sampleId, value.get().getSampleId());
    }
}
