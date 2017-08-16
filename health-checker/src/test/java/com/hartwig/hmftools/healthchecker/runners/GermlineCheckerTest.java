package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.HartwigException;
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
    private static final String SINGLE_SAMPLE_RUN = BASE_DIRECTORY + File.separator + "single_sample_run";

    private static final String SINGLE_SAMPLE = "sample";
    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final int EXPECTED_NUM_CHECKS_PER_SAMPLE = 5;

    private static final int EXPECTED_SINGLE_SAMPLE_SNPS = 2;
    private static final int EXPECTED_SINGLE_SAMPLE_INDELS = 0;
    private static final int EXPECTED_SINGLE_SAMPLE_HETEROZYGOUS_COUNT = 1;
    private static final int EXPECTED_SINGLE_SAMPLE_HETEROZYGOUS_COUNT_ABOVE_50 = 1;
    private static final int EXPECTED_SINGLE_SAMPLE_HETEROZYGOUS_COUNT_BELOW_50 = 0;

    private static final int EXPECTED_REF_SNPS = 55;
    private static final int EXPECTED_REF_INDELS = 4;
    private static final int EXPECTED_REF_HETEROZYGOUS_COUNT = 29;
    private static final int EXPECTED_REF_HETEROZYGOUS_COUNT_ABOVE_50 = 13;
    private static final int EXPECTED_REF_HETEROZYGOUS_COUNT_BELOW_50 = 16;

    private static final int EXPECTED_TUMOR_SNPS = 74;
    private static final int EXPECTED_TUMOR_INDELS = 4;
    private static final int EXPECTED_TUMOR_HETEROZYGOUS_COUNT = 62;
    private static final int EXPECTED_TUMOR_HETEROZYGOUS_COUNT_ABOVE_50 = 25;
    private static final int EXPECTED_TUMOR_HETEROZYGOUS_COUNT_BELOW_50 = 36;

    private final GermlineChecker checker = new GermlineChecker();

    @Test
    public void canCountSNPAndIndelsInSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(SINGLE_SAMPLE_RUN, SINGLE_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);

        assertEquals(CheckType.GERMLINE, result.getCheckType());
        final List<HealthCheck> checks = ((MultiValueResult) result).getChecks();
        assertChecks(checks, EXPECTED_SINGLE_SAMPLE_SNPS, EXPECTED_SINGLE_SAMPLE_INDELS, EXPECTED_SINGLE_SAMPLE_HETEROZYGOUS_COUNT,
                EXPECTED_SINGLE_SAMPLE_HETEROZYGOUS_COUNT_ABOVE_50, EXPECTED_SINGLE_SAMPLE_HETEROZYGOUS_COUNT_BELOW_50);
    }

    @Test
    public void canCountSNPAndIndelsInSomatic() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final PatientResult result = (PatientResult) checker.tryRun(runContext);

        assertEquals(CheckType.GERMLINE, result.getCheckType());
        final List<HealthCheck> refChecks = result.getRefSampleChecks();
        final List<HealthCheck> tumorChecks = result.getTumorSampleChecks();

        assertChecks(refChecks, EXPECTED_REF_SNPS, EXPECTED_REF_INDELS, EXPECTED_REF_HETEROZYGOUS_COUNT,
                EXPECTED_REF_HETEROZYGOUS_COUNT_ABOVE_50, EXPECTED_REF_HETEROZYGOUS_COUNT_BELOW_50);
        assertChecks(tumorChecks, EXPECTED_TUMOR_SNPS, EXPECTED_TUMOR_INDELS, EXPECTED_TUMOR_HETEROZYGOUS_COUNT,
                EXPECTED_TUMOR_HETEROZYGOUS_COUNT_ABOVE_50, EXPECTED_TUMOR_HETEROZYGOUS_COUNT_BELOW_50);
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

    private static void assertChecks(@NotNull final List<HealthCheck> checks, final int expectedCountSNP, final int expectedCountIndels,
            final int expectedHeterozygousCount, final int expectedHeterozygousCountAbove50VAF,
            final int expectedHeterozygousCountBelow50VAF) {
        assertEquals(EXPECTED_NUM_CHECKS_PER_SAMPLE, checks.size());

        assertCheck(checks, GermlineCheck.GERMLINE_SNP_COUNT, expectedCountSNP);
        assertCheck(checks, GermlineCheck.GERMLINE_INDEL_COUNT, expectedCountIndels);
        assertCheck(checks, GermlineCheck.GERMLINE_HETEROZYGOUS_COUNT, expectedHeterozygousCount);
        assertCheck(checks, GermlineCheck.GERMLINE_HETEROZYGOUS_COUNT_ABOVE_50_VAF, expectedHeterozygousCountAbove50VAF);
        assertCheck(checks, GermlineCheck.GERMLINE_HETEROZYGOUS_COUNT_BELOW_50_VAF, expectedHeterozygousCountBelow50VAF);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final GermlineCheck checkName,
            final int expectedValue) {
        final Optional<HealthCheck> optCheck =
                checks.stream().filter(check -> check.getCheckName().equals(checkName.toString())).findFirst();
        assert optCheck.isPresent();
        assertEquals(Integer.toString(expectedValue), optCheck.get().getValue());
    }
}
