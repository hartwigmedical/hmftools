package com.hartwig.hmftools.common.metrics;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricsFile {

    private static final String MEAN_COVERAGE_COLUMN = "MEAN_COVERAGE";
    private static final String COVERAGE_10X_COLUMN = "PCT_10X";
    private static final String COVERAGE_20X_COLUMN = "PCT_20X";
    private static final String COVERAGE_30X_COLUMN = "PCT_30X";
    private static final String COVERAGE_60X_COLUMN = "PCT_60X";

    private WGSMetricsFile() {
    }

    @NotNull
    public static WGSMetrics read(@NotNull String refFilePath) throws IOException {
        ImmutableWGSMetrics.Builder builder = ImmutableWGSMetrics.builder();

        appendRefValues(builder, refFilePath);

        return builder.build();
    }

    @NotNull
    public static WGSMetrics read(@NotNull String refFilePath, @NotNull String tumorFilePath) throws IOException {
        ImmutableWGSMetrics.Builder builder = ImmutableWGSMetrics.builder();

        appendRefValues(builder, refFilePath);
        appendTumorValues(builder, tumorFilePath);

        return builder.build();
    }

    private static void appendRefValues(@NotNull ImmutableWGSMetrics.Builder builder, @NotNull String refFilePath)
            throws IOException {
        WGSMetricsLines lines = WGSMetricsLines.fromFile(refFilePath);

        builder.refMeanCoverage(Double.parseDouble(lines.findValueByHeader(MEAN_COVERAGE_COLUMN)));
        builder.ref10xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_10X_COLUMN)));
        builder.ref20xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_20X_COLUMN)));
    }

    private static void appendTumorValues(@NotNull ImmutableWGSMetrics.Builder builder,
            @NotNull String tumorFilePath) throws IOException {
        WGSMetricsLines lines = WGSMetricsLines.fromFile(tumorFilePath);

        builder.tumorMeanCoverage(Double.parseDouble(lines.findValueByHeader(MEAN_COVERAGE_COLUMN)));
        builder.tumor30xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_30X_COLUMN)));
        builder.tumor60xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_60X_COLUMN)));
    }
}
