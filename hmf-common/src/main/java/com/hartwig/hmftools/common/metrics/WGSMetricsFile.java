package com.hartwig.hmftools.common.metrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.io.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.exception.MalformedFileException;
import com.hartwig.hmftools.common.io.path.PathPrefixSuffixFinder;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricsFile {

    // KODU: metrics files stores in {run}/QCStats/{sample}_dedup/{sample}_dedup_WGSMetrics.txt
    private static final String METRICS_BASE_DIRECTORY = "QCStats";
    private static final String METRICS_SUB_DIRECTORY_SUFFIX = "_dedup";
    private static final String METRICS_EXTENSION = "_WGSMetrics.txt";
    private static final String VALUE_SEPARATOR = "\t";

    private static final int MEAN_COVERAGE_INDEX = 1;
    private static final int COVERAGE_10X_INDEX = 13;
    private static final int COVERAGE_20X_INDEX = 15;
    private static final int COVERAGE_30X_INDEX = 17;
    private static final int COVERAGE_60X_INDEX = 20;

    private WGSMetricsFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull final String runDir, @NotNull final String sample) throws FileNotFoundException {
        String path = runDir + File.separator + METRICS_BASE_DIRECTORY + File.separator + sample + METRICS_SUB_DIRECTORY_SUFFIX;
        return PathPrefixSuffixFinder.build().findPath(path, sample, METRICS_EXTENSION).toString();
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
        String[] values = findValuesLine(refFilePath).split(VALUE_SEPARATOR);

        builder.refMeanCoverage(Double.parseDouble(values[MEAN_COVERAGE_INDEX]));
        builder.ref10xCoveragePercentage(Double.parseDouble(values[COVERAGE_10X_INDEX]));
        builder.ref20xCoveragePercentage(Double.parseDouble(values[COVERAGE_20X_INDEX]));

        return builder;
    }

    @NotNull
    private static ImmutableWGSMetrics.Builder appendTumorValues(@NotNull ImmutableWGSMetrics.Builder builder,
            @NotNull String tumorFilePath) throws IOException {
        String[] values = findValuesLine(tumorFilePath).split(VALUE_SEPARATOR);

        builder.tumorMeanCoverage(Double.parseDouble(values[MEAN_COVERAGE_INDEX]));
        builder.tumor30xCoveragePercentage(Double.parseDouble(values[COVERAGE_30X_INDEX]));
        builder.tumor60xCoveragePercentage(Double.parseDouble(values[COVERAGE_60X_INDEX]));

        return builder;
    }

    @NotNull
    private static String findValuesLine(@NotNull String filename) throws IOException {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if (lines.isEmpty()) {
            throw new EmptyFileException(filename);
        }
        final int index = findHeaderLineIndex(lines);
        if (index >= lines.size()) {
            throw new MalformedFileException(String.format("No value line found after header line in WGS Metrics file %s.", filename));
        }
        return lines.get(index + 1);
    }

    private static int findHeaderLineIndex(@NotNull final List<String> lines) throws MalformedFileException {
        final Optional<Integer> lineNumbers =
                IntStream.range(0, lines.size()).filter(index -> lines.get(index).contains("MEAN_COVERAGE")).boxed().findFirst();
        if (!lineNumbers.isPresent()) {
            throw new MalformedFileException(String.format("Could not find header line in WGS Metrics file with %s lines.", lines.size()));
        }
        return lineNumbers.get();
    }
}
