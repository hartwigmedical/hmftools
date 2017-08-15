package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.SomaticCheck;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class SomaticCheckerTest {

    private static final double EPSILON = 1.0e-4;
    private static final int EXPECTED_NUM_CHECKS = 9;

    private static final String RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("variants") + File.separator + "run";
    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final String MINIMAL_RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("variants") + File.separator + "run2";
    private static final String MINIMAL_REF_SAMPLE = "sample3";
    private static final String MINIMAL_TUMOR_SAMPLE = "sample4";

    private static final String AF_RUN_DIRECTORY = RunnerTestFunctions.getRunnerResourcePath("variants") + File.separator + "af_test";
    private static final String AF_REF_SAMPLE = "sample1";
    private static final String AF_TUMOR_SAMPLE = "sample2";

    private static final String INDELS = VariantType.INDEL.toString();
    private static final String SNP = VariantType.SNP.toString();

    private final SomaticChecker checker = new SomaticChecker();

    @Test
    public void canAnalyseTypicalMeltedVCF() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);

        final BaseResult result = checker.tryRun(runContext);
        final List<HealthCheck> checks = ((MultiValueResult) result).getChecks();

        Assert.assertEquals(CheckType.SOMATIC, result.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, checks.size());

        assertCheck(checks, SomaticCheck.COUNT_TOTAL.checkName(SNP), 987);
        assertCheck(checks, SomaticCheck.DBSNP_COUNT.checkName(SNP), 819);

        assertCheck(checks, SomaticCheck.COUNT_TOTAL.checkName(INDELS), 67);
        assertCheck(checks, SomaticCheck.DBSNP_COUNT.checkName(INDELS), 42);

        assertCheck(checks, SomaticCheck.PROPORTION_CHECK.checkName(SNP), 0.93643);
        assertCheck(checks, SomaticCheck.PROPORTION_CHECK.checkName(INDELS), 0.06356);
    }

    @Test
    public void runsCorrectlyForSingleSample() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final BaseResult result = checker.tryRun(runContext);
        assertTrue(result instanceof NoResult);
    }

    @Test
    public void errorYieldsCorrectOutputForSomatic() {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.errorRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
    }

    @Test
    public void errorYieldsCorrectOutputForSingleSample() {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final BaseResult result = checker.errorRun(runContext);
        assertTrue(result instanceof NoResult);
    }

    @Test
    public void canAnalyseMinimalVCF() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(MINIMAL_RUN_DIRECTORY, MINIMAL_REF_SAMPLE, MINIMAL_TUMOR_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.tryRun(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.getChecks().size());
    }

    @Test
    public void canDetermineAFs() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(AF_RUN_DIRECTORY, AF_REF_SAMPLE, AF_TUMOR_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.tryRun(runContext);
        final List<HealthCheck> checks = result.getChecks();
        assertCheck(checks, SomaticCheck.AF_LOWER_SD.checkName(), 0.1);
        assertCheck(checks, SomaticCheck.AF_MEDIAN.checkName(), 0.5);
        assertCheck(checks, SomaticCheck.AF_UPPER_SD.checkName(), 0.9);
    }

    @NotNull
    private static List<String> checkNames(@NotNull final List<HealthCheck> checks) {
        List<String> checkNames = Lists.newArrayList();
        for (HealthCheck check : checks) {
            checkNames.add(check.getCheckName());
        }
        return checkNames;
    }

    @Test
    public void errorGeneratesSameCheckNames() throws IOException, HartwigException {
        final RunContext normalRun = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);
        final List<HealthCheck> normalChecks = ((MultiValueResult) checker.tryRun(normalRun)).getChecks();
        final List<HealthCheck> errorChecks = ((MultiValueResult) checker.errorRun(normalRun)).getChecks();

        final List<String> normalCheckNames = checkNames(normalChecks);
        final List<String> errorCheckNames = checkNames(errorChecks);

        assertEquals(normalCheckNames, errorCheckNames);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String checkName, final double expectedValue) {
        final Optional<HealthCheck> report = checks.stream().filter(data -> data.getCheckName().equals(checkName)).findFirst();

        assert report.isPresent();
        final String check = report.get().getValue();
        double checkValue = Double.valueOf(check);
        assertEquals(expectedValue, checkValue, EPSILON);
    }
}
