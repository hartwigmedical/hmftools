package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class SomaticVariantsCheckerTest {

    private static final double EPSILON = 1.0e-4;
    private static final int EXPECTED_NUM_CHECKS = 4;

    private static final String BASE_DIRECTORY = Resources.getResource("somatics").getPath();
    private static final String RUN_DIRECTORY = BASE_DIRECTORY + File.separator + "run";
    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final String MINIMAL_RUN_DIRECTORY = BASE_DIRECTORY + File.separator + "run2";
    private static final String MINIMAL_REF_SAMPLE = "sample3";
    private static final String MINIMAL_TUMOR_SAMPLE = "sample4";

    private final SomaticVariantsChecker checker = new SomaticVariantsChecker();

    @Test
    public void canAnalyseTypicalSomaticVariantVCF() throws IOException, HartwigException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);

        final BaseResult result = checker.tryRun(runContext);
        final List<HealthCheck> checks = ((MultiValueResult) result).getChecks();

        Assert.assertEquals(CheckType.SOMATIC_VARIANTS, result.getCheckType());
        assertEquals(EXPECTED_NUM_CHECKS, checks.size());

        assertCheck(checks, SomaticVariantCheck.SOMATIC_SNP_COUNT.toString(), 991);
        assertCheck(checks, SomaticVariantCheck.SOMATIC_SNP_DBSNP_COUNT.toString(), 820);

        assertCheck(checks, SomaticVariantCheck.SOMATIC_INDEL_COUNT.toString(), 67);
        assertCheck(checks, SomaticVariantCheck.SOMATIC_INDEL_DBSNP_COUNT.toString(), 42);
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

    @NotNull
    private static List<String> checkNames(@NotNull final List<HealthCheck> checks) {
        final List<String> checkNames = Lists.newArrayList();
        for (final HealthCheck check : checks) {
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
