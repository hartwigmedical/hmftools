package com.hartwig.hmftools.healthchecker.runners;

import static com.google.common.collect.Iterables.getLast;

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
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
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
            "[EEE MMM d HH:mm:ss z yyyy]" + "[EEE MMM ppd HH:mm:ss z yyyy]" + "[EEE d MMM HH:mm:ss z yyyy]"
                    + "[EEE d MMM yyyy HH:mm:ss z]" + "[EEE MMM d yyyy HH:mm:ss z]";
    private static final String SOMATIC_LINE_TO_GET_DATE_FROM_REGEX = "End\\s+(Kinship|Finalize)";
    private static final String SINGLE_SAMPLE_LINE_TO_GET_DATE_FROM_REGEX = "End\\s+(germline variant annotation|Finalize)";
    private static final String DATE_LINE_FIELD_SEPARATOR = "\t";

    private static final String PIPELINE_LOG_REGEX = "PipelineCheck.log";
    private static final String PIPELINE_VERSION_REGEX = "Pipeline version:";
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

        return toFinalResult(runContext, runDate, pipelineVersion);
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
        return toFinalResult(runContext, HealthCheckConstants.ERROR_VALUE, pipelineVersion);
    }

    @NotNull
    private BaseResult toFinalResult(@NotNull final RunContext runContext, @NotNull final String runDate,
            @NotNull final String pipelineVersion) {
        final List<HealthCheck> refMetaData = toHealthCheckList(runContext.refSample(), runContext, runDate,
                pipelineVersion);
        HealthCheck.log(LOGGER, refMetaData);

        if (runContext.isSomaticRun()) {
            final List<HealthCheck> tumorMetaData = toHealthCheckList(runContext.tumorSample(), runContext, runDate,
                    pipelineVersion);

            HealthCheck.log(LOGGER, tumorMetaData);

            return new PatientResult(checkType(), refMetaData, tumorMetaData);
        } else {
            return new MultiValueResult(checkType(), refMetaData);
        }
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
        final Predicate<String> dateLineFilter = runContext.isSomaticRun() ?
                doesLineMatch(SOMATIC_LINE_TO_GET_DATE_FROM_REGEX) :
                doesLineMatch(SINGLE_SAMPLE_LINE_TO_GET_DATE_FROM_REGEX);
        List<String> dateLines = LineReader.build().readLines(dateTimeLogPath, dateLineFilter);
        final String date = datePart(getLast(dateLines).split(DATE_LINE_FIELD_SEPARATOR));
        final DateTimeFormatter inFormatter = DateTimeFormatter.ofPattern(DATE_IN_FORMAT, Locale.ENGLISH);
        final LocalDateTime formattedDate = LocalDateTime.parse(date, inFormatter);
        final DateTimeFormatter outFormatter = DateTimeFormatter.ofPattern(DATE_OUT_FORMAT, Locale.ENGLISH);
        return outFormatter.format(formattedDate);
    }

    @NotNull
    private static String datePart(@NotNull String[] parts) {
        if (parts[0].contains(" ")) {
            return parts[1].trim();
        } else {
            return parts[2].trim();
        }
    }

    @NotNull
    private static String extractPipelineVersion(@NotNull final String runDirectory)
            throws IOException, HartwigException {
        final Path pipelineLogPath = PathRegexFinder.build().findPath(runDirectory, PIPELINE_LOG_REGEX);
        final List<String> versionsLine = LineReader.build().readLines(pipelineLogPath,
                doesLineMatch(PIPELINE_VERSION_REGEX));
        return versionsLine.get(PIPELINE_VERSION_LINE_INDEX).split(PIPELINE_VERSION_LINE_SEPARATOR)[1].trim();
    }

    @NotNull
    private static Predicate<String> doesLineMatch(@NotNull final String regex) {
        return new Predicate<String>() {
            @Override
            public boolean test(final String line) {
                return line.matches(String.format("%s.*", regex));
            }

            @Override
            public String toString() {
                return String.format("doesLineMatch %s", regex);
            }
        };
    }
}
