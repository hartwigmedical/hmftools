package com.hartwig.hmftools.healthchecker.runners;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.runners.checks.CopynumberCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.exception.EmptyFileException;
import com.hartwig.hmftools.healthchecker.exception.HealthChecksException;
import com.hartwig.hmftools.healthchecker.exception.MalformedFileException;
import com.hartwig.hmftools.healthchecker.io.dir.RunContext;
import com.hartwig.hmftools.healthchecker.io.dir.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class CopynumberCheckerTest {
    private static final String RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("copynumber");

    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";
    private static final String MALFORMED_SAMPLE = "sample3";
    private static final String NO_CNVS_FOUND_SAMPLE = "sample4";
    private static final String EMPTY_RATIOS_SAMPLE = "sample5";
    private static final String NO_RATIOS_SAMPLE = "sample6";
    private static final String MISSING_CNVS_SAMPLE = "sample7";
    private static final int EXPECTED_NUM_CHECKS = 2;
    private static final long EXPECTED_GAIN = 252;
    private static final long EXPECTED_LOSS = 11561;

    private final CopynumberChecker checker = new CopynumberChecker();

    @Test
    public void correctInputYieldsCorrectOutput() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);
        assertResult(result, EXPECTED_GAIN, EXPECTED_LOSS);
    }

    @Test
    public void errorYieldsCorrectOutput() {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
    }

    @Test(expected = MalformedFileException.class)
    public void noGainLossTagsYieldMalformedException() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, REF_SAMPLE, MALFORMED_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test
    public void emptyCopyNumberVariantsAllowedIfRatiosPresent() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, NO_CNVS_FOUND_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);
        assertResult(result, 0, 0);
    }

    @Test(expected = EmptyFileException.class)
    public void emptyCopyNumberVariantsNotAllowedIfRatiosEmpty() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, EMPTY_RATIOS_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = FileNotFoundException.class)
    public void emptyCopyNumberVariantsNotAllowedIfRatiosMissing() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, NO_RATIOS_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    @Test(expected = FileNotFoundException.class)
    public void missingCopyNumberVariantsNotAllowedEvenIfRatiosPresent() throws IOException, HealthChecksException {
        final RunContext runContext = TestRunContextFactory.forTest(RUN_DIRECTORY, MISSING_CNVS_SAMPLE, TUMOR_SAMPLE);
        checker.tryRun(runContext);
    }

    private static void assertResult(@NotNull final BaseResult baseResult, long expectedGain, long expectedLoss) {
        final MultiValueResult result = (MultiValueResult) baseResult;
        Assert.assertEquals(CheckType.COPYNUMBER, result.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());

        final HealthCheck gainCheck = extractHealthCheck(result.getChecks(), CopynumberCheck.COPYNUMBER_GENOME_GAIN);

        assertEquals(TUMOR_SAMPLE, gainCheck.getSampleId());
        assertEquals(Long.toString(expectedGain), gainCheck.getValue());

        final HealthCheck lossCheck = extractHealthCheck(result.getChecks(), CopynumberCheck.COPYNUMBER_GENOME_LOSS);

        assertEquals(TUMOR_SAMPLE, lossCheck.getSampleId());
        assertEquals(Long.toString(expectedLoss), lossCheck.getValue());
    }

    @NotNull
    private static HealthCheck extractHealthCheck(@NotNull final List<HealthCheck> checks,
            @NotNull final CopynumberCheck checkName) {
        final Optional<HealthCheck> optCheck = checks.stream().filter(
                check -> check.getCheckName().equals(checkName.toString())).findFirst();

        assert optCheck.isPresent();
        return optCheck.get();
    }
}