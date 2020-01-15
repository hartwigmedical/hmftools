package com.hartwig.hmftools.common.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.healthchecker.result.QCValue;
import com.hartwig.hmftools.common.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class MetricsChecker implements HealthChecker {

    @NotNull
    private final String refSample;
    @Nullable
    private final String tumorSample;
    @NotNull
    private final String metricsDirectory;

    public MetricsChecker(@NotNull final String refSample, @Nullable final String tumorSample, @NotNull final String metricsDirectory) {
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
    private WGSMetrics extractMetrics() throws IOException {
        String refFile = WGSMetricsFile.generateFilename(metricsDirectory, refSample);
        if (tumorSample != null) {
            String tumorFile = WGSMetricsFile.generateFilename(metricsDirectory, tumorSample);
            return WGSMetricsFile.read(refFile, tumorFile);
        } else {
            return WGSMetricsFile.read(refFile);
        }
    }
}
