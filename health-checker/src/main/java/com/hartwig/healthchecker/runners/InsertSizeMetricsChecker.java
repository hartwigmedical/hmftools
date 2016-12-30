package com.hartwig.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.healthchecker.runners.checks.HealthCheck;
import com.hartwig.healthchecker.runners.checks.InsertSizeMetricsCheck;
import com.hartwig.healthchecker.exception.HealthChecksException;
import com.hartwig.healthchecker.exception.LineNotFoundException;
import com.hartwig.healthchecker.io.dir.RunContext;
import com.hartwig.healthchecker.io.path.PathPrefixSuffixFinder;
import com.hartwig.healthchecker.io.reader.FileReader;
import com.hartwig.healthchecker.resource.ResourceWrapper;
import com.hartwig.healthchecker.result.BaseResult;
import com.hartwig.healthchecker.result.PatientResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.INSERT_SIZE)
public class InsertSizeMetricsChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(InsertSizeMetricsChecker.class);

    // KODU: metrics files stores in {run}/QCStats/{sample}_dedup/{sample}<>.insert_size_metrics
    private static final String METRICS_BASE_DIRECTORY = "QCStats";
    private static final String METRICS_SUB_DIRECTORY_SUFFIX = "_dedup";
    private static final String INSERT_SIZE_METRICS_EXTENSION = ".insert_size_metrics";
    private static final String VALUE_SEPARATOR = "\t";

    public InsertSizeMetricsChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.INSERT_SIZE;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HealthChecksException {
        final List<HealthCheck> refChecks = extractChecksForSample(runContext.runDirectory(), runContext.refSample());
        final List<HealthCheck> tumorChecks = extractChecksForSample(runContext.runDirectory(),
                runContext.tumorSample());

        return toPatientResult(refChecks, tumorChecks);
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        return toPatientResult(getErrorChecksForSample(runContext.refSample()),
                getErrorChecksForSample(runContext.tumorSample()));
    }

    @NotNull
    private static List<HealthCheck> getErrorChecksForSample(@NotNull final String sampleId) {
        final List<HealthCheck> checks = Lists.newArrayList();
        for (InsertSizeMetricsCheck check : InsertSizeMetricsCheck.values()) {
            checks.add(new HealthCheck(sampleId, check.toString(), HealthCheckConstants.ERROR_VALUE));
        }
        return checks;
    }

    @NotNull
    private BaseResult toPatientResult(@NotNull final List<HealthCheck> refChecks,
            @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(checkType(), refChecks, tumorChecks);
    }

    @NotNull
    private static List<HealthCheck> extractChecksForSample(@NotNull final String runDirectory,
            @NotNull final String sampleId) throws IOException, HealthChecksException {
        final String basePath = getBasePathForSample(runDirectory, sampleId);
        final Path insertSizeMetricsPath = PathPrefixSuffixFinder.build().findPath(basePath, sampleId,
                INSERT_SIZE_METRICS_EXTENSION);
        final List<String> lines = FileReader.build().readLines(insertSizeMetricsPath);

        final HealthCheck medianReport = getCheck(insertSizeMetricsPath.toString(), lines, sampleId,
                InsertSizeMetricsCheck.MAPPING_MEDIAN_INSERT_SIZE);
        final HealthCheck width70PerReport = getCheck(insertSizeMetricsPath.toString(), lines, sampleId,
                InsertSizeMetricsCheck.MAPPING_WIDTH_OF_70_PERCENT);
        return Arrays.asList(medianReport, width70PerReport);
    }

    @NotNull
    private static String getBasePathForSample(@NotNull final String runDirectory, @NotNull final String sampleId) {
        return runDirectory + File.separator + METRICS_BASE_DIRECTORY + File.separator + sampleId
                + METRICS_SUB_DIRECTORY_SUFFIX;
    }

    @NotNull
    private static HealthCheck getCheck(@NotNull final String filePath, @NotNull final List<String> lines,
            @NotNull final String sampleId, @NotNull final InsertSizeMetricsCheck check) throws LineNotFoundException {
        final String value = getValueFromLine(filePath, lines, check.getFieldName(), check.getColumnIndex());
        return new HealthCheck(sampleId, check.toString(), value);
    }

    @NotNull
    private static String getValueFromLine(@NotNull final String filePath, @NotNull final List<String> lines,
            @NotNull final String filter, final int fieldIndex) throws LineNotFoundException {
        final int index = findLineIndex(filePath, lines, filter);
        final String line = lines.get(index + 1);
        final String[] lineValues = line.split(VALUE_SEPARATOR);
        return lineValues[fieldIndex];
    }

    private static int findLineIndex(@NotNull final String filePath, @NotNull final List<String> lines,
            @NotNull final String filter) throws LineNotFoundException {
        final Optional<Integer> lineNumbers = IntStream.range(0, lines.size()).filter(
                index -> lines.get(index).contains(filter)).mapToObj(index -> index).findFirst();
        if (!lineNumbers.isPresent()) {
            throw new LineNotFoundException(filePath, filter);
        }
        return lineNumbers.get();
    }
}
