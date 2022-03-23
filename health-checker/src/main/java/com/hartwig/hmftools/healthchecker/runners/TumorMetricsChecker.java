package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;

public class TumorMetricsChecker extends FileBasedHealthChecker<WGSMetrics> {

    public TumorMetricsChecker(@NotNull String metricsFile) throws IOException {
        super(WGSMetricsFile.read(metricsFile),
                wgsMetrics -> List.of(ImmutableQCValue.builder()
                                .type(QCValueType.TUM_COVERAGE_30X)
                                .value(String.valueOf(wgsMetrics.coverage30xPercentage()))
                                .build(),
                        ImmutableQCValue.builder()
                                .type(QCValueType.TUM_COVERAGE_60X)
                                .value(String.valueOf(wgsMetrics.coverage60xPercentage()))
                                .build()));
    }
}
