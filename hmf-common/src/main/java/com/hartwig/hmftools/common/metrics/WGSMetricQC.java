package com.hartwig.hmftools.common.metrics;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricQC {

    public static final double MIN_REF_10X_COVERAGE = 0.9;
    public static final double MIN_REF_20X_COVERAGE = 0.7;
    public static final double MIN_TUMOR_30X_COVERAGE = 0.8;
    public static final double MIN_TUMOR_60X_COVERAGE = 0.65;

    private WGSMetricQC() {
    }

    @NotNull
    public static ImmutableWGSMetricWithQC buildWithQCMetric(@NotNull WGSMetrics metrics) {
        // This function is only expected to be called for somatic runs.
        Double tumor30xCoveragePercentage = metrics.tumor30xCoveragePercentage();
        Double tumor60xCoveragePercentage = metrics.tumor60xCoveragePercentage();

        assert tumor30xCoveragePercentage != null;
        assert tumor60xCoveragePercentage != null;

        boolean wgsQCRef10 = metrics.ref10xCoveragePercentage() >= MIN_REF_10X_COVERAGE;
        boolean wgsQCRef20 = metrics.ref20xCoveragePercentage() >= MIN_REF_20X_COVERAGE;
        boolean wgsQCTumor30 = tumor30xCoveragePercentage >= MIN_TUMOR_30X_COVERAGE;
        boolean wgsQCTumor60 = tumor60xCoveragePercentage >= MIN_TUMOR_60X_COVERAGE;

        boolean wgsQC = wgsQCRef10 && wgsQCRef20 && wgsQCTumor30 && wgsQCTumor60;

        return ImmutableWGSMetricWithQC.builder().wgsMetrics(metrics).qcMetric(wgsQC).build();
    }
}
