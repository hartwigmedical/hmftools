package com.hartwig.hmftools.common.metrics;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricQC {

    private WGSMetricQC() {
    }

    public static ImmutableWGSMetricWithQC checkQCMetric(@NotNull WGSMetrics metrics) {

        double ref10xCoverage = metrics.ref10xCoveragePercentage();
        double ref20xCoverage = metrics.ref20xCoveragePercentage();

        // We only write metrics for somatic runs.
        assert metrics.tumor30xCoveragePercentage() != null;
        assert metrics.tumor60xCoveragePercentage() != null;

        Double tumor30xCoverage = metrics.tumor30xCoveragePercentage();
        Double tumor60xCoverage = metrics.tumor60xCoveragePercentage();

        boolean WGSqc = true;

        if (tumor30xCoverage != null && tumor60xCoverage != null) {
            if (ref10xCoverage < 0.9 || ref20xCoverage < 0.7 || tumor30xCoverage < 0.8 || tumor60xCoverage < 0.65) {
                WGSqc = false;
            }
        }

        return ImmutableWGSMetricWithQC.builder().wgsMetrics(metrics).qcMetric(WGSqc).build();
    }
}
