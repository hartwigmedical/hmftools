package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.context.TestRunContextFactory;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class StrelkaCheckerTest {

    private static final double EPSILON = 1.0e-4;
    private static final int EXPECTED_NUM_CHECKS = 6;

    private static final String BASE_DIRECTORY = Resources.getResource("strelka").getPath();
    private static final String RUN_DIRECTORY = BASE_DIRECTORY + File.separator + "run";
    private static final String REF_SAMPLE = "sample1";
    private static final String TUMOR_SAMPLE = "sample2";

    private static final String MINIMAL_RUN_DIRECTORY = BASE_DIRECTORY + File.separator + "run2";
    private static final String MINIMAL_REF_SAMPLE = "sample3";
    private static final String MINIMAL_TUMOR_SAMPLE = "sample4";

    private final StrelkaChecker checker = new StrelkaChecker();

    @Test
    public void canAnalyseTypicalSomaticVariantVCF() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(RUN_DIRECTORY, REF_SAMPLE, TUMOR_SAMPLE);

        final BaseResult result = checker.run(runContext);
        final List<HealthCheck> checks = ((MultiValueResult) result).checks();

        Assert.assertEquals(CheckType.STRELKA, result.checkType());
        assertEquals(EXPECTED_NUM_CHECKS, checks.size());

        assertCheck(checks, StrelkaCheck.SOMATIC_SNP_COUNT.toString(), 990);
        assertCheck(checks, StrelkaCheck.SOMATIC_SNP_DBSNP_COUNT.toString(), 820);

        assertCheck(checks, StrelkaCheck.SOMATIC_INDEL_COUNT.toString(), 67);
        assertCheck(checks, StrelkaCheck.SOMATIC_INDEL_DBSNP_COUNT.toString(), 42);

        assertCheck(checks, StrelkaCheck.SOMATIC_MNP_COUNT.toString(), 1);
        assertCheck(checks, StrelkaCheck.SOMATIC_MNP_DBSNP_COUNT.toString(), 0);
    }

    @Test
    public void runsCorrectlyForSingleSample() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSingleSampleTest(RUN_DIRECTORY, REF_SAMPLE);
        final BaseResult result = checker.run(runContext);
        assertTrue(result instanceof NoResult);
    }

    @Test
    public void canAnalyseMinimalVCF() throws IOException {
        final RunContext runContext = TestRunContextFactory.forSomaticTest(MINIMAL_RUN_DIRECTORY, MINIMAL_REF_SAMPLE, MINIMAL_TUMOR_SAMPLE);
        final MultiValueResult result = (MultiValueResult) checker.run(runContext);
        assertEquals(EXPECTED_NUM_CHECKS, result.checks().size());
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String checkName, final double expectedValue) {
        final Optional<HealthCheck> report = checks.stream().filter(data -> data.getCheckName().equals(checkName)).findFirst();

        assert report.isPresent();
        final String check = report.get().getValue();
        double checkValue = Double.valueOf(check);
        assertEquals(expectedValue, checkValue, EPSILON);
    }
}
