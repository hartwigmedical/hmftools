package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class MetricsChecker implements HealthChecker {

    @NotNull
    private final String refWgsMetricsFile;
    @Nullable
    private final String tumWgsMetricsFile;

    public MetricsChecker(@NotNull final String refWgsMetricsFile, @Nullable final String tumWgsMetricsFile) {
        this.refWgsMetricsFile = refWgsMetricsFile;
        this.tumWgsMetricsFile = tumWgsMetricsFile;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {
        WGSMetrics metrics = extractMetrics();

        List<QCValue> qcValues = Lists.newArrayList();
        qcValues.add(ImmutableQCValue.builder()
                .type(QCValueType.REF_COVERAGE_10X)
                .value(String.valueOf(metrics.ref10xCoveragePercentage()))
                .build());
        qcValues.add(ImmutableQCValue.builder()
                .type(QCValueType.REF_COVERAGE_20X)
                .value(String.valueOf(metrics.ref20xCoveragePercentage()))
                .build());

        if (tumWgsMetricsFile != null) {
            assert metrics.tumor30xCoveragePercentage() != null;
            assert metrics.tumor60xCoveragePercentage() != null;

            qcValues.add(ImmutableQCValue.builder()
                    .type(QCValueType.TUM_COVERAGE_30X)
                    .value(String.valueOf(metrics.tumor30xCoveragePercentage()))
                    .build());
            qcValues.add(ImmutableQCValue.builder()
                    .type(QCValueType.TUM_COVERAGE_60X)
                    .value(String.valueOf(metrics.tumor60xCoveragePercentage()))
                    .build());
        }

        return qcValues;
    }

    @NotNull
    private WGSMetrics extractMetrics() throws IOException {
        if (tumWgsMetricsFile != null) {
            return WGSMetricsFile.read(refWgsMetricsFile, tumWgsMetricsFile);
        } else {
            return WGSMetricsFile.read(refWgsMetricsFile);
        }
    }
}
