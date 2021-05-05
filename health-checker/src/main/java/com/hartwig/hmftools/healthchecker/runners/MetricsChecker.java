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
    private final String tumorWgsMetricsFile;

    public MetricsChecker(@NotNull final String refWgsMetricsFile, @Nullable final String tumorWgsMetricsFile) {
        this.refWgsMetricsFile = refWgsMetricsFile;
        this.tumorWgsMetricsFile = tumorWgsMetricsFile;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {
        WGSMetrics refMetrics = WGSMetricsFile.read(refWgsMetricsFile);

        List<QCValue> qcValues = Lists.newArrayList();
        qcValues.add(ImmutableQCValue.builder()
                .type(QCValueType.REF_COVERAGE_10X)
                .value(String.valueOf(refMetrics.coverage10xPercentage()))
                .build());
        qcValues.add(ImmutableQCValue.builder()
                .type(QCValueType.REF_COVERAGE_20X)
                .value(String.valueOf(refMetrics.coverage20xPercentage()))
                .build());

        if (tumorWgsMetricsFile != null) {
            WGSMetrics tumorMetrics = WGSMetricsFile.read(tumorWgsMetricsFile);

            qcValues.add(ImmutableQCValue.builder()
                    .type(QCValueType.TUM_COVERAGE_30X)
                    .value(String.valueOf(tumorMetrics.coverage30xPercentage()))
                    .build());
            qcValues.add(ImmutableQCValue.builder()
                    .type(QCValueType.TUM_COVERAGE_60X)
                    .value(String.valueOf(tumorMetrics.coverage60xPercentage()))
                    .build());
        }

        return qcValues;
    }
}
