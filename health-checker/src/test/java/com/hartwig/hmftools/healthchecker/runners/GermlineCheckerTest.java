package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.context.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.GermlineCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineCheckerTest {

    private static final String BASE_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("variants");
    private static final String RUN_DIRECTORY = BASE_DIRECTORY + File.separator + "run";
    private static final String RUN_DIRECTORY_V1_10 = BASE_DIRECTORY + File.separator + "run_v1_10";
    private static final String RUN_DIRECTORY_V1_9 = BASE_DIRECTORY + File.separator + "run_v1_9";
    private static final String SINGLE_SAMPLE_RUN = BASE_DIRECTORY + File.separator + "single_sample_run";

    private static final String SINGLE_SAMPLE = "sample";
    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final int EXPECTED_NUM_CHECKS_PER_SAMPLE = 2;
    private static final int EXPECTED_SINGLE_SAMPLE_SNPS = 2;
    private static final int EXPECTED_SINGLE_SAMPLE_INDELS = 0;
    private static final int EXPECTED_REF_SNPS = 55;
    private static final int EXPECTED_REF_INDELS = 4;
    private static final int EXPECTED_TUMOR_SNPS = 74;
    private static final int EXPECTED_TUMOR_INDELS = 4;

    private final GermlineChecker checker = new GermlineChecker();

    @Test
    public void canCountSNPAndIndelsInSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(SINGLE_SAMPLE_RUN, SINGLE_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        assertEquals(CheckType.GERMLINE, result.getCheckType());
        final List<HealthCheck> checks = ((MultiValueResult) result).getChecks();
        assertChecks(checks, EXPECTED_SINGLE_SAMPLE_SNPS, EXPECTED_SINGLE_SAMPLE_INDELS);
    }

    @Test
    public void canCountSNPAndIndelsInSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        assertEquals(CheckType.GERMLINE, result.getCheckType());
        final List<HealthCheck> refChecks = ((PatientResult) result).getRefSampleChecks();
        final List<HealthCheck> tumorChecks = ((PatientResult) result).getTumorSampleChecks();

        assertChecks(refChecks, EXPECTED_REF_SNPS, EXPECTED_REF_INDELS);
        assertChecks(tumorChecks, EXPECTED_TUMOR_SNPS, EXPECTED_TUMOR_INDELS);
    }

    @Test
    public void canCountSNPAndIndels_V1_10() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY_V1_10, REF_SAMPLE,
                TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.tryRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, result.getTumorSampleChecks().size());
        assertChecks(result.getRefSampleChecks(), EXPECTED_REF_SNPS, EXPECTED_REF_INDELS);
        assertChecks(result.getTumorSampleChecks(), EXPECTED_TUMOR_SNPS, EXPECTED_TUMOR_INDELS);
    }

    @Test
    public void canCountSNPAndIndelsV1_9InSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY_V1_9, REF_SAMPLE,
                TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.tryRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, result.getTumorSampleChecks().size());
        assertChecks(result.getRefSampleChecks(), 1, 0);
        assertChecks(result.getTumorSampleChecks(), 1, 0);
    }

    @Test
    public void errorYieldsCorrectOutputForSomatic() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, result.getRefSampleChecks().size());
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, result.getTumorSampleChecks().size());
    }

    @Test
    public void errorYieldsCorrectOutputForSingleSample() {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, result.getChecks().size());
    }

    private static void assertChecks(@NotNull final List<HealthCheck> checks, final long expectedCountSNP,
            final long expectedCountIndels) {
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, checks.size());

        final Optional<HealthCheck> snpCheck = checks.stream().filter(
                check -> check.getCheckName().equals(GermlineCheck.VARIANTS_GERMLINE_SNP.toString())).findFirst();
        assert snpCheck.isPresent();
        assertEquals(Long.toString(expectedCountSNP), snpCheck.get().getValue());

        final Optional<HealthCheck> indelCheck = checks.stream().filter(
                check -> check.getCheckName().equals(GermlineCheck.VARIANTS_GERMLINE_INDELS.toString())).findFirst();
        assert indelCheck.isPresent();
        assertEquals(Long.toString(expectedCountIndels), indelCheck.get().getValue());
    }
}
