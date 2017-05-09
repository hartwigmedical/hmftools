package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.exception.LineNotFoundException;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.SummaryMetricsCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SummaryMetricsCheckerTest {

    private static final String RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("metrics");

    private static final int EXPECTED_NUM_CHECKS = 5;

    private static final String REF_SAMPLE = "sample1";
    private static final String REF_PF_MISMATCH_RATE = "0.006024";
    private static final String REF_PF_INDEL_RATE = "0.000261";
    private static final String REF_STRAND_BALANCE = "0.399972";
    private static final String REF_PCT_CHIMERA = "0.000212";
    private static final String REF_PCT_ADAPTER = "0.000046";

    private static final String TUMOR_SAMPLE = "sample2";
    private static final String TUMOR_PF_MISMATCH_RATE = "0.005024";
    private static final String TUMOR_PF_INDEL_RATE = "0.000262";
    private static final String TUMOR_STRAND_BALANCE = "0.499972";
    private static final String TUMOR_PCT_CHIMERA = "0.000112";
    private static final String TUMOR_PCT_ADAPTER = "0.000056";

    private static final String EMPTY_SAMPLE = "sample3";
    private static final String INCORRECT_SAMPLE = "sample4";
    private static final String NON_EXISTING_SAMPLE = "sample5";

    private final SummaryMetricsChecker checker = new SummaryMetricsChecker();

    @Test
    public void correctInputYieldsCorrectOutputForSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);
        assertSomaticResult(result);
    }

    @Test
    public void correctInputYieldsCorrectOutputForSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.tryRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
        assertCheck(result.getChecks(), REF_SAMPLE, SummaryMetricsCheck.MAPPING_PCT_ADAPTER, REF_PCT_ADAPTER);
        assertCheck(result.getChecks(), REF_SAMPLE, SummaryMetricsCheck.MAPPING_PCT_CHIMERA, REF_PCT_CHIMERA);
        assertCheck(result.getChecks(), REF_SAMPLE, SummaryMetricsCheck.MAPPING_PF_INDEL_RATE, REF_PF_INDEL_RATE);
        assertCheck(result.getChecks(), REF_SAMPLE, SummaryMetricsCheck.MAPPING_PF_MISMATCH_RATE,
                REF_PF_MISMATCH_RATE);
        assertCheck(result.getChecks(), REF_SAMPLE, SummaryMetricsCheck.MAPPING_STRAND_BALANCE, REF_STRAND_BALANCE);
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

        assertEquals(CheckType.SUMMARY_METRICS, patientResult.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS, patientResult.getTumorSampleChecks().size());

        assertSomaticField(patientResult, SummaryMetricsCheck.MAPPING_PF_MISMATCH_RATE, REF_PF_MISMATCH_RATE,
                TUMOR_PF_MISMATCH_RATE);
        assertSomaticField(patientResult, SummaryMetricsCheck.MAPPING_PF_INDEL_RATE, REF_PF_INDEL_RATE,
                TUMOR_PF_INDEL_RATE);
        assertSomaticField(patientResult, SummaryMetricsCheck.MAPPING_STRAND_BALANCE, REF_STRAND_BALANCE,
                TUMOR_STRAND_BALANCE);
        assertSomaticField(patientResult, SummaryMetricsCheck.MAPPING_PCT_CHIMERA, REF_PCT_CHIMERA, TUMOR_PCT_CHIMERA);
        assertSomaticField(patientResult, SummaryMetricsCheck.MAPPING_PCT_ADAPTER, REF_PCT_ADAPTER, TUMOR_PCT_ADAPTER);
    }

    private static void assertSomaticField(@NotNull final PatientResult result,
            @NotNull final SummaryMetricsCheck field, @NotNull final String refValue, @NotNull final String tumValue) {
        assertCheck(result.getRefSampleChecks(), REF_SAMPLE, field, refValue);
        assertCheck(result.getTumorSampleChecks(), TUMOR_SAMPLE, field, tumValue);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String sampleId,
            @NotNull final SummaryMetricsCheck check, @NotNull final String expectedValue) {
        final Optional<HealthCheck> value = checks.stream().filter(
                p -> p.getCheckName().equals(check.toString())).findFirst();
        assert value.isPresent();

        assertEquals(expectedValue, value.get().getValue());
        assertEquals(sampleId, value.get().getSampleId());
    }
}
