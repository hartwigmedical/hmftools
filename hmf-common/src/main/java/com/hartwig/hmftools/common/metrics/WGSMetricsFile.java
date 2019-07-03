package com.hartwig.hmftools.common.metrics;

import java.io.FileNotFoundException;
import java.io.IOException;

import com.hartwig.hmftools.common.io.path.PathPrefixSuffixFinder;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricsFile {

    private static final String METRICS_EXTENSION = ".wgsmetrics";

    private static final String MEAN_COVERAGE_COLUMN = "MEAN_COVERAGE";
    private static final String COVERAGE_10X_COLUMN = "PCT_10X";
    private static final String COVERAGE_20X_COLUMN = "PCT_20X";
    private static final String COVERAGE_30X_COLUMN = "PCT_30X";
    private static final String COVERAGE_60X_COLUMN = "PCT_60X";

    private WGSMetricsFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull final String metricsDir, @NotNull final String sample) throws FileNotFoundException {
        return PathPrefixSuffixFinder.build().findPath(metricsDir, sample, METRICS_EXTENSION).toString();
    }

    @NotNull
    public static WGSMetrics read(@NotNull String refFilePath) throws IOException {
        ImmutableWGSMetrics.Builder builder = ImmutableWGSMetrics.builder();

        builder = appendRefValues(builder, refFilePath);

        return builder.build();
    }

    @NotNull
    public static WGSMetrics read(@NotNull String refFilePath, @NotNull String tumorFilePath) throws IOException {
        ImmutableWGSMetrics.Builder builder = ImmutableWGSMetrics.builder();

        builder = appendRefValues(builder, refFilePath);
        builder = appendTumorValues(builder, tumorFilePath);

        return builder.build();
    }

    @NotNull
    private static ImmutableWGSMetrics.Builder appendRefValues(@NotNull ImmutableWGSMetrics.Builder builder, @NotNull String refFilePath)
            throws IOException {
        WGSMetricsLines lines = WGSMetricsLines.fromFile(refFilePath);

        builder.refMeanCoverage(Double.parseDouble(lines.findValueByHeader(MEAN_COVERAGE_COLUMN)));
        builder.ref10xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_10X_COLUMN)));
        builder.ref20xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_20X_COLUMN)));

        return builder;
    }

    @NotNull
    private static ImmutableWGSMetrics.Builder appendTumorValues(@NotNull ImmutableWGSMetrics.Builder builder,
            @NotNull String tumorFilePath) throws IOException {
        WGSMetricsLines lines = WGSMetricsLines.fromFile(tumorFilePath);

        builder.tumorMeanCoverage(Double.parseDouble(lines.findValueByHeader(MEAN_COVERAGE_COLUMN)));
        builder.tumor30xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_30X_COLUMN)));
        builder.tumor60xCoveragePercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_60X_COLUMN)));

        return builder;
    }
}
