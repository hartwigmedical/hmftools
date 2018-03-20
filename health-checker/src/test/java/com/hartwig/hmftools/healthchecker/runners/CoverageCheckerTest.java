package com.hartwig.hmftools.healthchecker.runners;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.metrics.ImmutableWGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CoverageCheckerTest {

    private static final String REF_SAMPLE = "sample1";
    private static final int REF_NUMBER_OF_CHECKS = 2;
    private static final double REF_COVERAGE_10X = 0.1;
    private static final double REF_COVERAGE_20X = 0.2;

    private static final String TUMOR_SAMPLE = "sample2";
    private static final int TUMOR_NUMBER_OF_CHECKS = 2;
    private static final double TUMOR_COVERAGE_30X = 0.3;
    private static final double TUMOR_COVERAGE_60X = 0.4;

    @Test
    public void worksForSomaticRun() {
        WGSMetrics metrics = ImmutableWGSMetrics.builder()
                .ref10xCoveragePercentage(REF_COVERAGE_10X)
                .ref20xCoveragePercentage(REF_COVERAGE_20X)
                .refMeanCoverage(0D)
                .tumor30xCoveragePercentage(TUMOR_COVERAGE_30X)
                .tumor60xCoveragePercentage(TUMOR_COVERAGE_60X)
                .build();

        PatientResult result = (PatientResult) CoverageChecker.toCheckResult(metrics, REF_SAMPLE, TUMOR_SAMPLE);

        List<HealthCheck> refResult = result.refSampleChecks();
        assertEquals(REF_NUMBER_OF_CHECKS, refResult.size());
        assertCheck(refResult, REF_SAMPLE, CoverageCheck.COVERAGE_10X, REF_COVERAGE_10X);
        assertCheck(refResult, REF_SAMPLE, CoverageCheck.COVERAGE_20X, REF_COVERAGE_20X);

        List<HealthCheck> tumorResult = result.tumorSampleChecks();
        assertEquals(TUMOR_NUMBER_OF_CHECKS, tumorResult.size());
        assertCheck(tumorResult, TUMOR_SAMPLE, CoverageCheck.COVERAGE_30X, TUMOR_COVERAGE_30X);
        assertCheck(tumorResult, TUMOR_SAMPLE, CoverageCheck.COVERAGE_60X, TUMOR_COVERAGE_60X);
    }

    @Test
    public void worksForSingleSampleRun() {
        WGSMetrics metrics = ImmutableWGSMetrics.builder()
                .ref10xCoveragePercentage(REF_COVERAGE_10X)
                .ref20xCoveragePercentage(REF_COVERAGE_20X)
                .refMeanCoverage(0D)
                .build();

        MultiValueResult result = (MultiValueResult) CoverageChecker.toCheckResult(metrics, REF_SAMPLE, null);

        List<HealthCheck> refResult = result.checks();
        assertEquals(REF_NUMBER_OF_CHECKS, refResult.size());
        assertCheck(refResult, REF_SAMPLE, CoverageCheck.COVERAGE_10X, REF_COVERAGE_10X);
        assertCheck(refResult, REF_SAMPLE, CoverageCheck.COVERAGE_20X, REF_COVERAGE_20X);
    }

    private static void assertCheck(@NotNull final List<HealthCheck> checks, @NotNull final String sampleId,
            @NotNull final CoverageCheck check, final double expectedValue) {
        final Optional<HealthCheck> value = checks.stream().filter(p -> p.getCheckName().equals(check.toString())).findFirst();
        assert value.isPresent();

        assertEquals(expectedValue, Double.parseDouble(value.get().getValue()), 1e-10);
        assertEquals(sampleId, value.get().getSampleId());
    }
}
