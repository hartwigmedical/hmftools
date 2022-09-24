package com.hartwig.hmftools.common.metrics;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricQC
{
    public static final double MIN_REF_10X_COVERAGE = 0.9;
    public static final double MIN_REF_20X_COVERAGE = 0.7;
    public static final double MIN_TUMOR_30X_COVERAGE = 0.8;
    public static final double MIN_TUMOR_60X_COVERAGE = 0.65;

    @NotNull
    public static WGSMetricWithQC buildWithQCMetric(final WGSMetrics refMetrics, final WGSMetrics tumorMetrics)
    {
        boolean wgsQCRef10 = refMetrics.coverage10xPercentage() >= MIN_REF_10X_COVERAGE;
        boolean wgsQCRef20 = refMetrics.coverage20xPercentage() >= MIN_REF_20X_COVERAGE;
        boolean wgsQCTumor30 = tumorMetrics.coverage30xPercentage() >= MIN_TUMOR_30X_COVERAGE;
        boolean wgsQCTumor60 = tumorMetrics.coverage60xPercentage() >= MIN_TUMOR_60X_COVERAGE;

        boolean sufficientCoverage = wgsQCRef10 && wgsQCRef20 && wgsQCTumor30 && wgsQCTumor60;

        return ImmutableWGSMetricWithQC.builder()
                .refMetrics(refMetrics)
                .tumorMetrics(tumorMetrics)
                .sufficientCoverage(sufficientCoverage)
                .build();
    }
}
