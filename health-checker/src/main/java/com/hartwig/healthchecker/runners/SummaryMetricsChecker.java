package com.hartwig.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.healthchecker.runners.checks.HealthCheck;
import com.hartwig.healthchecker.runners.checks.SummaryMetricsCheck;
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
@ResourceWrapper(type = CheckType.SUMMARY_METRICS)
public class SummaryMetricsChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(SummaryMetricsChecker.class);

    // KODU: metrics files stores in {run}/QCStats/{sample}_dedup/{sample}<>.alignment_summary_metrics
    private static final String METRICS_BASE_DIRECTORY = "QCStats";
    private static final String METRICS_SUB_DIRECTORY_SUFFIX = "_dedup";
    private static final String ALIGNMENT_SUMMARY_METRICS_EXTENSION = ".alignment_summary_metrics";

    private static final String PICARD_CATEGORY_TO_READ = "PAIR";
    private static final String VALUE_SEPARATOR = "\t";

    public SummaryMetricsChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.SUMMARY_METRICS;
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
        for (SummaryMetricsCheck check : SummaryMetricsCheck.values()) {
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
        Path alignmentSummaryMetricsPath = PathPrefixSuffixFinder.build().findPath(basePath, sampleId,
                ALIGNMENT_SUMMARY_METRICS_EXTENSION);
        final List<String> lines = FileReader.build().readLines(alignmentSummaryMetricsPath);

        final Optional<String> searchedLine = lines.stream().filter(
                fileLine -> fileLine.startsWith(PICARD_CATEGORY_TO_READ)).findFirst();
        if (!searchedLine.isPresent()) {
            throw new LineNotFoundException(alignmentSummaryMetricsPath.toString(), PICARD_CATEGORY_TO_READ);
        }

        final HealthCheck pfIndelRate = getCheck(searchedLine.get(), sampleId,
                SummaryMetricsCheck.MAPPING_PF_INDEL_RATE);
        final HealthCheck pctAdapter = getCheck(searchedLine.get(), sampleId, SummaryMetricsCheck.MAPPING_PCT_ADAPTER);
        final HealthCheck pctChimeras = getCheck(searchedLine.get(), sampleId,
                SummaryMetricsCheck.MAPPING_PCT_CHIMERA);
        final HealthCheck pfMisMatch = getCheck(searchedLine.get(), sampleId,
                SummaryMetricsCheck.MAPPING_PF_MISMATCH_RATE);
        final HealthCheck strandBalance = getCheck(searchedLine.get(), sampleId,
                SummaryMetricsCheck.MAPPING_STRAND_BALANCE);
        return Arrays.asList(pfIndelRate, pctAdapter, pctChimeras, pfMisMatch, strandBalance);
    }

    @NotNull
    private static String getBasePathForSample(@NotNull final String runDirectory, @NotNull final String sampleId) {
        return runDirectory + File.separator + METRICS_BASE_DIRECTORY + File.separator + sampleId
                + METRICS_SUB_DIRECTORY_SUFFIX;
    }

    @NotNull
    private static HealthCheck getCheck(@NotNull final String line, @NotNull final String sampleId,
            @NotNull final SummaryMetricsCheck check) {
        final String value = line.split(VALUE_SEPARATOR)[check.getIndex()];
        return new HealthCheck(sampleId, check.toString(), value);
    }
}
