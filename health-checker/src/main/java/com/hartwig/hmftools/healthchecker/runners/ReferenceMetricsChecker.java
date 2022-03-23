package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;

public class ReferenceMetricsChecker extends FileBasedHealthChecker<WGSMetrics> {

    public ReferenceMetricsChecker(@NotNull String metricsFile) throws IOException {
        super(WGSMetricsFile.read(metricsFile),
                wgsMetrics -> List.of(ImmutableQCValue.builder()
                                .type(QCValueType.REF_COVERAGE_10X)
                                .value(String.valueOf(wgsMetrics.coverage10xPercentage()))
                                .build(),
                        ImmutableQCValue.builder()
                                .type(QCValueType.REF_COVERAGE_20X)
                                .value(String.valueOf(wgsMetrics.coverage20xPercentage()))
                                .build()));
    }
}
