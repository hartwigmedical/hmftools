package com.hartwig.hmftools.common.metrics;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricsTestFactory {

    private WGSMetricsTestFactory() {
    }

    @NotNull
    public static WGSMetrics createMinimalTestWGSMetrics() {
        return ImmutableWGSMetrics.builder()
                .meanCoverage(0D)
                .sdCoverage(0D)
                .medianCoverage(0)
                .madCoverage(0)
                .pctExcMapQ(0D)
                .pctExcDupe(0D)
                .pctExcUnpaired(0D)
                .pctExcBaseQ(0D)
                .pctExcOverlap(0D)
                .pctExcCapped(0D)
                .pctExcTotal(0D)
                .coverage10xPercentage(0D)
                .coverage20xPercentage(0D)
                .coverage30xPercentage(0D)
                .coverage60xPercentage(0D)
                .build();
    }
}
