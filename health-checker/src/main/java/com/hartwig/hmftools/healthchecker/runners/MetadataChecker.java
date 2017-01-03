package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathRegexFinder;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.MetadataCheck;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.METADATA)
public class MetadataChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(MetadataChecker.class);

    private static final String LOG_FILENAME_FORMAT = "%s.log";
    private static final String DATE_OUT_FORMAT = "yyyy-MM-dd";
    private static final String DATE_IN_FORMAT =
            "[EEE MMM d HH:mm:ss z yyyy]" + "[EEE d MMM HH:mm:ss z yyyy]" + "[EEE d MMM yyyy HH:mm:ss z]"
                    + "[EEE MMM d yyyy HH:mm:ss z]";
    private static final String REGEX_SPLIT = "\t";
    private static final String LINE_TO_GET_DATE_FROM = "End Kinship";
    private static final int DATE_LINE_INDEX = 0;

    private static final String PIPELINE_LOG_REGEX = "PipelineCheck.log";
    private static final String PIPELINE_VERSION = "Pipeline version:";
    private static final String PIPELINE_VERSION_LINE_SEPARATOR = ":";
    private static final int PIPELINE_VERSION_LINE_INDEX = 0;

    public MetadataChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.METADATA;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final String runDate = extractRunDate(runContext);
        final String pipelineVersion = extractPipelineVersion(runContext.runDirectory());

        return toPatientResult(runContext, runDate, pipelineVersion);
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        String pipelineVersion = HealthCheckConstants.ERROR_VALUE;
        try {
            pipelineVersion = extractPipelineVersion(runContext.runDirectory());
        } catch (IOException | HartwigException exception) {
            // KODU: This is to work around the fact that pipeline version and run date are linked.
            // If run date extraction fails, pipeline version wont be available.
            // Should be two separate checkers...
        }
        return toPatientResult(runContext, HealthCheckConstants.ERROR_VALUE, pipelineVersion);
    }

    @NotNull
    private BaseResult toPatientResult(@NotNull final RunContext runContext, @NotNull final String runDate,
            @NotNull final String pipelineVersion) {
        final List<HealthCheck> refMetaData = toHealthCheckList(runContext.refSample(), runContext, runDate,
                pipelineVersion);
        final List<HealthCheck> tumorMetaData = toHealthCheckList(runContext.tumorSample(), runContext, runDate,
                pipelineVersion);

        HealthCheck.log(LOGGER, refMetaData);
        HealthCheck.log(LOGGER, tumorMetaData);

        return new PatientResult(checkType(), refMetaData, tumorMetaData);
    }

    @NotNull
    private static List<HealthCheck> toHealthCheckList(@NotNull final String sampleId,
            @NotNull final RunContext runContext, @NotNull final String runDate,
            @NotNull final String pipelineVersion) {
        return Lists.newArrayList(new HealthCheck(sampleId, MetadataCheck.RUN_NAME.toString(), runContext.runName()),
                new HealthCheck(sampleId, MetadataCheck.HAS_PASSED_TESTS.toString(),
                        Boolean.toString(runContext.hasPassedTests())),
                new HealthCheck(sampleId, MetadataCheck.RUN_DATE.toString(), runDate),
                new HealthCheck(sampleId, MetadataCheck.PIPELINE_VERSION.toString(), pipelineVersion));
    }

    @NotNull
    private static String extractRunDate(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final Path dateTimeLogPath = PathRegexFinder.build().findPath(runContext.runDirectory(),
                String.format(LOG_FILENAME_FORMAT, runContext.runName()));
        final List<String> dateLine = LineReader.build().readLines(dateTimeLogPath,
                doesLineStartWith(LINE_TO_GET_DATE_FROM));
        // KODU: Replacing all double spaces with single spaces to solve issue mentioned in run2 of MetadataCheckerTest
        final String date = dateLine.get(DATE_LINE_INDEX).split(REGEX_SPLIT)[1].trim().replaceAll("  ", " ");
        final DateTimeFormatter inFormatter = DateTimeFormatter.ofPattern(DATE_IN_FORMAT, Locale.ENGLISH);
        final LocalDateTime formattedDate = LocalDateTime.parse(date, inFormatter);
        final DateTimeFormatter outFormatter = DateTimeFormatter.ofPattern(DATE_OUT_FORMAT, Locale.ENGLISH);
        return outFormatter.format(formattedDate);
    }

    @NotNull
    private static String extractPipelineVersion(@NotNull final String runDirectory)
            throws IOException, HartwigException {
        final Path pipelineLogPath = PathRegexFinder.build().findPath(runDirectory, PIPELINE_LOG_REGEX);
        final List<String> versionsLine = LineReader.build().readLines(pipelineLogPath,
                doesLineStartWith(PIPELINE_VERSION));
        return versionsLine.get(PIPELINE_VERSION_LINE_INDEX).split(PIPELINE_VERSION_LINE_SEPARATOR)[1].trim();
    }

    @NotNull
    private static Predicate<String> doesLineStartWith(@NotNull final String prefix) {
        return line -> line.startsWith(prefix);
    }
}
