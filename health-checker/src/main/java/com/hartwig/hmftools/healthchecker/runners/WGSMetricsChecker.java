package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.common.exception.LineNotFoundException;
import com.hartwig.hmftools.common.io.dir.RunContext;
import com.hartwig.hmftools.common.io.path.PathPrefixSuffixFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.WGSMetricsCheck;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.WGS_METRICS)
public class WGSMetricsChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(WGSMetricsChecker.class);

    // KODU: metrics files stores in {run}/QCStats/{sample}_dedup/{sample}_dedup_WGSMetrics.txt
    private static final String METRICS_BASE_DIRECTORY = "QCStats";
    private static final String METRICS_SUB_DIRECTORY_SUFFIX = "_dedup";
    private static final String WGS_METRICS_EXTENSION = "_WGSMetrics.txt";
    private static final String VALUE_SEPARATOR = "\t";

    public WGSMetricsChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.WGS_METRICS;
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
        for (WGSMetricsCheck check : WGSMetricsCheck.values()) {
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
        final Path wgsMetricsPath = PathPrefixSuffixFinder.build().findPath(basePath, sampleId, WGS_METRICS_EXTENSION);
        final List<String> lines = FileReader.build().readLines(wgsMetricsPath);

        final HealthCheck coverageMean = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_MEAN);
        final HealthCheck coverageMedian = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_MEDIAN);
        final HealthCheck coverageSD = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_SD);
        final HealthCheck coverage10X = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_10X);
        final HealthCheck coverage20X = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_20X);
        final HealthCheck coverage30X = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_30X);
        final HealthCheck coverage60X = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_60X);
        final HealthCheck coverageBaseQ = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_PCT_EXC_BASEQ);
        final HealthCheck coverageDupe = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_PCT_EXC_DUPE);
        final HealthCheck coverageMapQ = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_PCT_EXC_MAPQ);
        final HealthCheck coverageOverlap = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_PCT_EXC_OVERLAP);
        final HealthCheck coverageTotal = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_PCT_EXC_TOTAL);
        final HealthCheck coverageUnpaired = getCheck(wgsMetricsPath.toString(), lines, sampleId,
                WGSMetricsCheck.COVERAGE_PCT_EXC_UNPAIRED);

        return Arrays.asList(coverageMean, coverageMedian, coverageSD, coverage10X, coverage20X, coverage30X,
                coverage60X, coverageBaseQ, coverageDupe, coverageMapQ, coverageOverlap, coverageTotal,
                coverageUnpaired);
    }

    @NotNull
    private static String getBasePathForSample(@NotNull final String runDirectory, @NotNull final String sampleId) {
        return runDirectory + File.separator + METRICS_BASE_DIRECTORY + File.separator + sampleId
                + METRICS_SUB_DIRECTORY_SUFFIX;
    }

    @NotNull
    private static HealthCheck getCheck(@NotNull final String filePath, @NotNull final List<String> lines,
            @NotNull final String sampleId, @NotNull final WGSMetricsCheck check) throws LineNotFoundException {
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
