package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.exception.LineNotFoundException;
import com.hartwig.hmftools.common.io.path.PathPrefixSuffixFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CoverageChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(CoverageChecker.class);

    // KODU: metrics files stores in {run}/QCStats/{sample}_dedup/{sample}_dedup_WGSMetrics.txt
    private static final String METRICS_BASE_DIRECTORY = "QCStats";
    private static final String METRICS_SUB_DIRECTORY_SUFFIX = "_dedup";
    private static final String WGS_METRICS_EXTENSION = "_WGSMetrics.txt";
    private static final String VALUE_SEPARATOR = "\t";

    public CoverageChecker() {
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final List<HealthCheck> refChecks = extractChecksForSample(runContext.runDirectory(), runContext.refSample());
        if (runContext.isSomaticRun()) {
            final List<HealthCheck> tumorChecks = extractChecksForSample(runContext.runDirectory(), runContext.tumorSample());

            return toPatientResult(refChecks, tumorChecks);
        } else {
            return toMultiValueResult(refChecks);
        }
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        if (runContext.isSomaticRun()) {
            return toPatientResult(getErrorChecksForSample(runContext.refSample()), getErrorChecksForSample(runContext.tumorSample()));
        } else {
            return toMultiValueResult(getErrorChecksForSample(runContext.refSample()));
        }
    }

    @NotNull
    private static List<HealthCheck> getErrorChecksForSample(@NotNull final String sampleId) {
        final List<HealthCheck> checks = Lists.newArrayList();
        for (final CoverageCheck check : CoverageCheck.values()) {
            checks.add(new HealthCheck(sampleId, check.toString(), HealthCheckConstants.ERROR_VALUE));
        }
        return checks;
    }

    @NotNull
    private BaseResult toPatientResult(@NotNull final List<HealthCheck> refChecks, @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(CheckType.COVERAGE, refChecks, tumorChecks);
    }

    @NotNull
    private static BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);

        return new MultiValueResult(CheckType.COVERAGE, checks);
    }

    @NotNull
    private static List<HealthCheck> extractChecksForSample(@NotNull final String runDirectory, @NotNull final String sampleId)
            throws IOException, HartwigException {
        final String basePath = getBasePathForSample(runDirectory, sampleId);
        final Path wgsMetricsPath = PathPrefixSuffixFinder.build().findPath(basePath, sampleId, WGS_METRICS_EXTENSION);
        final List<String> lines = FileReader.build().readLines(wgsMetricsPath);

        final HealthCheck coverage10X = getCheck(wgsMetricsPath.toString(), lines, sampleId, CoverageCheck.COVERAGE_10X);
        final HealthCheck coverage20X = getCheck(wgsMetricsPath.toString(), lines, sampleId, CoverageCheck.COVERAGE_20X);
        final HealthCheck coverage30X = getCheck(wgsMetricsPath.toString(), lines, sampleId, CoverageCheck.COVERAGE_30X);
        final HealthCheck coverage60X = getCheck(wgsMetricsPath.toString(), lines, sampleId, CoverageCheck.COVERAGE_60X);

        return Arrays.asList(coverage10X, coverage20X, coverage30X, coverage60X);
    }

    @NotNull
    private static String getBasePathForSample(@NotNull final String runDirectory, @NotNull final String sampleId) {
        return runDirectory + File.separator + METRICS_BASE_DIRECTORY + File.separator + sampleId + METRICS_SUB_DIRECTORY_SUFFIX;
    }

    @NotNull
    private static HealthCheck getCheck(@NotNull final String filePath, @NotNull final List<String> lines, @NotNull final String sampleId,
            @NotNull final CoverageCheck check) throws LineNotFoundException {
        final String value = getValueFromLine(filePath, lines, check.getFieldName(), check.getColumnIndex());
        return new HealthCheck(sampleId, check.toString(), value);
    }

    @NotNull
    private static String getValueFromLine(@NotNull final String filePath, @NotNull final List<String> lines, @NotNull final String filter,
            final int fieldIndex) throws LineNotFoundException {
        final int index = findLineIndex(filePath, lines, filter);
        final String line = lines.get(index + 1);
        final String[] lineValues = line.split(VALUE_SEPARATOR);
        return lineValues[fieldIndex];
    }

    private static int findLineIndex(@NotNull final String filePath, @NotNull final List<String> lines, @NotNull final String filter)
            throws LineNotFoundException {
        final Optional<Integer> lineNumbers =
                IntStream.range(0, lines.size()).filter(index -> lines.get(index).contains(filter)).boxed().findFirst();
        if (!lineNumbers.isPresent()) {
            throw new LineNotFoundException(filePath, filter);
        }
        return lineNumbers.get();
    }
}
