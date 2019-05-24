package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CoverageChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(CoverageChecker.class);

    @NotNull
    private final String refSample;
    @Nullable
    private final String tumorSample;
    @NotNull
    private final String metricsDirectory;

    public CoverageChecker(@NotNull final String refSample, @Nullable final String tumorSample, @NotNull final String metricsDirectory) {
        this.refSample = refSample;
        this.tumorSample = tumorSample;
        this.metricsDirectory = metricsDirectory;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {
        final WGSMetrics metrics = extractMetrics();

        List<QCValue> qcValues = Lists.newArrayList();
        qcValues.add(ImmutableQCValue.of(QCValueType.REF_COVERAGE_10X, String.valueOf(metrics.ref10xCoveragePercentage())));
        qcValues.add(ImmutableQCValue.of(QCValueType.REF_COVERAGE_20X, String.valueOf(metrics.ref20xCoveragePercentage())));

        if (tumorSample != null) {
            assert metrics.tumor30xCoveragePercentage() != null;
            assert metrics.tumor60xCoveragePercentage() != null;

            qcValues.add(ImmutableQCValue.of(QCValueType.TUMOR_COVERAGE_30X, String.valueOf(metrics.tumor30xCoveragePercentage())));
            qcValues.add(ImmutableQCValue.of(QCValueType.TUMOR_COVERAGE_60X, String.valueOf(metrics.tumor60xCoveragePercentage())));
        }

        return qcValues;
    }

    @NotNull
    @VisibleForTesting
    static BaseResult toCheckResult(@NotNull WGSMetrics metrics, @NotNull String refSample, @Nullable String tumorSample) {
        final List<HealthCheck> refChecks = Lists.newArrayList(new HealthCheck(refSample,
                        CoverageCheck.COVERAGE_10X.toString(),
                        String.valueOf(metrics.ref10xCoveragePercentage())),
                new HealthCheck(refSample, CoverageCheck.COVERAGE_20X.toString(), String.valueOf(metrics.ref20xCoveragePercentage())));

        if (tumorSample != null) {
            assert metrics.tumor30xCoveragePercentage() != null;
            assert metrics.tumor60xCoveragePercentage() != null;

            final List<HealthCheck> tumorChecks = Lists.newArrayList(new HealthCheck(tumorSample,
                            CoverageCheck.COVERAGE_30X.toString(),
                            String.valueOf(metrics.tumor30xCoveragePercentage())),
                    new HealthCheck(tumorSample,
                            CoverageCheck.COVERAGE_60X.toString(),
                            String.valueOf(metrics.tumor60xCoveragePercentage())));
            return toPatientResult(refChecks, tumorChecks);
        } else {
            return toMultiValueResult(refChecks);
        }
    }

    @NotNull
    private WGSMetrics extractMetrics() throws IOException {
        String refFile = WGSMetricsFile.generateFilenamePv5(metricsDirectory, refSample);
        if (tumorSample != null) {
            String tumorFile = WGSMetricsFile.generateFilenamePv5(metricsDirectory, tumorSample);
            return WGSMetricsFile.read(refFile, tumorFile);
        } else {
            return WGSMetricsFile.read(refFile);
        }
    }

    @NotNull
    private static BaseResult toPatientResult(@NotNull final List<HealthCheck> refChecks, @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(CheckType.COVERAGE, refChecks, tumorChecks);
    }

    @NotNull
    private static BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);

        return new MultiValueResult(CheckType.COVERAGE, checks);
    }
}
