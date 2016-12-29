package com.hartwig.healthchecker.checks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.healthchecker.common.checks.CheckType;
import com.hartwig.healthchecker.common.checks.HealthCheck;
import com.hartwig.healthchecker.common.exception.EmptyFileException;
import com.hartwig.healthchecker.common.exception.HealthChecksException;
import com.hartwig.healthchecker.common.exception.LineNotFoundException;
import com.hartwig.healthchecker.common.exception.MalformedFileException;
import com.hartwig.healthchecker.common.io.dir.RunContext;
import com.hartwig.healthchecker.common.io.dir.TestRunContextFactory;
import com.hartwig.healthchecker.common.result.BaseResult;
import com.hartwig.healthchecker.common.result.PatientResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RealignerCheckerTest {

    private static final String REF_CHANGED_READS_PROPORTION = "0.04000";
    private static final String TUMOR_CHANGED_READS_PROPORTION = "0.04444";

    private static final String RUN_DIRECTORY = Resources.getResource("checks/realigner").getPath();

    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final String EMPTY_SAMPLE = "sample3";
    private static final String INCORRECT_SAMPLE = "sample4";
    private static final String NON_EXISTING_SAMPLE = "sample5";
    private static final String MALFORMED_SAMPLE = "sample6";

    private final RealignerChecker checker = new RealignerChecker();

    @Test
    public void correctInputYieldsCorrectOutput() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult report = checker.tryRun(runContext);
        assertReport(report);
    }

    @Test
    public void errorRunYieldsCorrectNumberOfChecks() {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(1, result.getRefSampleChecks().size());
        assertEquals(1, result.getTumorSampleChecks().size());
    }

    @Test(expected = EmptyFileException.class)
    public void emptyFileYieldsEmptyFileException() throws IOException, HealthChecksException {
        RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, EMPTY_SAMPLE, EMPTY_SAMPLE);
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

    @Test(expected = MalformedFileException.class)
    public void malformedLineYieldsMalformedFileException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, MALFORMED_SAMPLE, MALFORMED_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertReport(@NotNull final BaseResult result) {
        assertEquals(CheckType.REALIGNER, result.getCheckType());
        assertNotNull(result);
        assertField(result, RealignerChecker.REALIGNER_CHECK_NAME, REF_CHANGED_READS_PROPORTION,
                TUMOR_CHANGED_READS_PROPORTION);
    }

    private static void assertField(@NotNull final BaseResult result, @NotNull final String field,
            @NotNull final String refValue, @NotNull final String tumValue) {
        assertBaseData(((PatientResult) result).getRefSampleChecks(), REF_SAMPLE, field, refValue);
        assertBaseData(((PatientResult) result).getTumorSampleChecks(), TUMOR_SAMPLE, field, tumValue);
    }

    private static void assertBaseData(@NotNull final List<HealthCheck> checks, @NotNull final String sampleId,
            @NotNull final String checkName, @NotNull final String expectedValue) {
        final Optional<HealthCheck> value = checks.stream().filter(
                p -> p.getCheckName().equals(checkName)).findFirst();
        assert value.isPresent();

        assertEquals(expectedValue, value.get().getValue());
        assertEquals(sampleId, value.get().getSampleId());
    }
}
