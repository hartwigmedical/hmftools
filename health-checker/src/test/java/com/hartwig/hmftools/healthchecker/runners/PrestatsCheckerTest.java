package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.PrestatsCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PrestatsCheckerTest {

    private static final String RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("prestats");

    private static final int EXPECTED_CHECKS_NUM = 34;

    private static final String REF_SAMPLE = "sample1";
    private static final int REF_TOTAL_SEQUENCES = 700;

    private static final String TUMOR_SAMPLE = "sample2";
    private static final int TUMOR_TOTAL_SEQUENCES = 1450;

    private static final String EMPTY_TOTAL_SEQUENCE_SAMPLE = "sample3";
    private static final String INCOMPLETE_SAMPLE = "sample4";
    private static final String NON_EXISTING_SAMPLE = "sample5";

    @NotNull
    private final PrestatsChecker checker = new PrestatsChecker();

    @Test
    public void correctInputYieldsCorrectOutputForSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        assertEquals(CheckType.PRESTATS, result.getCheckType());
        assertRefChecks(((PatientResult) result).getRefSampleChecks());
        assertTumorChecks(((PatientResult) result).getTumorSampleChecks());
    }

    @Test
    public void correctInputYieldsCorrectOutputForSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        assertEquals(CheckType.PRESTATS, result.getCheckType());
        assertRefChecks(((MultiValueResult) result).getChecks());
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecksForSomatic() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_CHECKS_NUM, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_CHECKS_NUM, result.getTumorSampleChecks().size());
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecksForSingleSample() {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_CHECKS_NUM, result.getChecks().size());
    }

    @Test(expected = EmptyFileException.class)
    public void emptyTotalSequenceFileYieldsEmptyFileException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, EMPTY_TOTAL_SEQUENCE_SAMPLE,
                EMPTY_TOTAL_SEQUENCE_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test
    public void incompleteInputYieldsIncompleteOutput() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, INCOMPLETE_SAMPLE,
                INCOMPLETE_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);
        final List<HealthCheck> refResults = ((PatientResult) result).getRefSampleChecks();

        assertEquals(EXPECTED_CHECKS_NUM, refResults.size());
        assertResult(refResults, PrestatsCheck.PRESTATS_SEQUENCE_DUPLICATION_LEVELS, PrestatsCheckValue.FAIL, 0);
    }

    @Test(expected = IOException.class)
    public void nonExistingFileYieldsIOException() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, NON_EXISTING_SAMPLE,
                NON_EXISTING_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertRefChecks(@NotNull final List<HealthCheck> checks) {
        assertEquals(EXPECTED_CHECKS_NUM, checks.size());
        for (HealthCheck check : checks) {
            assertEquals(REF_SAMPLE, check.getSampleId());
        }
        assertResult(checks, PrestatsCheck.PRESTATS_NUMBER_OF_READS.toString(), REF_TOTAL_SEQUENCES);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_QUALITY, PrestatsCheckValue.FAIL, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_QUALITY, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_QUALITY, PrestatsCheckValue.PASS, 2);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_TILE_SEQUENCE_QUALITY, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_TILE_SEQUENCE_QUALITY, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_TILE_SEQUENCE_QUALITY, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_QUALITY_SCORES, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_QUALITY_SCORES, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_QUALITY_SCORES, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_CONTENT, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_CONTENT, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_CONTENT, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_GC_CONTENT, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_GC_CONTENT, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_GC_CONTENT, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_N_CONTENT, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_N_CONTENT, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_N_CONTENT, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_LENGTH_DISTRIBUTION, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_LENGTH_DISTRIBUTION, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_LENGTH_DISTRIBUTION, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_DUPLICATION_LEVELS, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_DUPLICATION_LEVELS, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_DUPLICATION_LEVELS, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_OVERREPRESENTED_SEQUENCES, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_OVERREPRESENTED_SEQUENCES, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_OVERREPRESENTED_SEQUENCES, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_ADAPTER_CONTENT, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_ADAPTER_CONTENT, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_ADAPTER_CONTENT, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_KMER_CONTENT, PrestatsCheckValue.FAIL, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_KMER_CONTENT, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_KMER_CONTENT, PrestatsCheckValue.PASS, 2);
    }

    private static void assertTumorChecks(@NotNull final List<HealthCheck> checks) {
        assertEquals(EXPECTED_CHECKS_NUM, checks.size());
        for (HealthCheck check : checks) {
            assertEquals(TUMOR_SAMPLE, check.getSampleId());
        }
        assertResult(checks, PrestatsCheck.PRESTATS_NUMBER_OF_READS.toString(), TUMOR_TOTAL_SEQUENCES);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_QUALITY, PrestatsCheckValue.FAIL, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_QUALITY, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_QUALITY, PrestatsCheckValue.PASS, 2);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_TILE_SEQUENCE_QUALITY, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_TILE_SEQUENCE_QUALITY, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_TILE_SEQUENCE_QUALITY, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_QUALITY_SCORES, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_QUALITY_SCORES, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_QUALITY_SCORES, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_CONTENT, PrestatsCheckValue.FAIL, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_CONTENT, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_SEQUENCE_CONTENT, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_GC_CONTENT, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_GC_CONTENT, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_SEQUENCE_GC_CONTENT, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_N_CONTENT, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_N_CONTENT, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_PER_BASE_N_CONTENT, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_LENGTH_DISTRIBUTION, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_LENGTH_DISTRIBUTION, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_LENGTH_DISTRIBUTION, PrestatsCheckValue.PASS, 4);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_DUPLICATION_LEVELS, PrestatsCheckValue.FAIL, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_DUPLICATION_LEVELS, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_SEQUENCE_DUPLICATION_LEVELS, PrestatsCheckValue.PASS, 3);
        assertResult(checks, PrestatsCheck.PRESTATS_OVERREPRESENTED_SEQUENCES, PrestatsCheckValue.FAIL, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_OVERREPRESENTED_SEQUENCES, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_OVERREPRESENTED_SEQUENCES, PrestatsCheckValue.PASS, 2);
        assertResult(checks, PrestatsCheck.PRESTATS_ADAPTER_CONTENT, PrestatsCheckValue.FAIL, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_ADAPTER_CONTENT, PrestatsCheckValue.WARN, 1);
        assertResult(checks, PrestatsCheck.PRESTATS_ADAPTER_CONTENT, PrestatsCheckValue.PASS, 2);
        assertResult(checks, PrestatsCheck.PRESTATS_KMER_CONTENT, PrestatsCheckValue.FAIL, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_KMER_CONTENT, PrestatsCheckValue.WARN, 0);
        assertResult(checks, PrestatsCheck.PRESTATS_KMER_CONTENT, PrestatsCheckValue.PASS, 4);
    }

    private static void assertResult(@NotNull final List<HealthCheck> checks, @NotNull final PrestatsCheck check,
            @NotNull final PrestatsCheckValue value, final int expectedCount) {
        final String checkName = PrestatsChecker.toCheckName(check, value);
        assertResult(checks, checkName, expectedCount);
    }

    private static void assertResult(@NotNull final List<HealthCheck> checks, @NotNull final String checkName,
            final int expectedValue) {
        Optional<HealthCheck> optCheck = checks.stream().filter(
                p -> p.getCheckName().equals(checkName)).findFirst();

        assert optCheck.isPresent();
        final String actualValue = optCheck.get().getValue();

        assertEquals(Integer.toString(expectedValue), actualValue);
    }
}
