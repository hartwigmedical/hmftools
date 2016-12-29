package com.hartwig.healthchecker.checks;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.healthchecker.common.checks.CheckType;
import com.hartwig.healthchecker.common.checks.ErrorHandlingChecker;
import com.hartwig.healthchecker.common.checks.HealthCheck;
import com.hartwig.healthchecker.common.checks.HealthCheckConstants;
import com.hartwig.healthchecker.common.checks.HealthChecker;
import com.hartwig.healthchecker.common.exception.HealthChecksException;
import com.hartwig.healthchecker.common.exception.LineNotFoundException;
import com.hartwig.healthchecker.common.exception.MalformedFileException;
import com.hartwig.healthchecker.common.io.dir.RunContext;
import com.hartwig.healthchecker.common.io.path.PathPrefixSuffixFinder;
import com.hartwig.healthchecker.common.io.reader.FileReader;
import com.hartwig.healthchecker.common.resource.ResourceWrapper;
import com.hartwig.healthchecker.common.result.BaseResult;
import com.hartwig.healthchecker.common.result.PatientResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.REALIGNER)
public class RealignerChecker extends ErrorHandlingChecker implements HealthChecker {

    @VisibleForTesting
    static final String REALIGNER_CHECK_NAME = "MAPPING_REALIGNER_CHANGED_ALIGNMENTS";

    private static final Logger LOGGER = LogManager.getLogger(RealignerChecker.class);

    private static final String REALIGNER_CHECK_PRECISION = "#0.00000";
    private static final String MALFORMED_FILE_MSG = "Malformed %s path was expecting %s in file";
    private static final String REALIGNER_BASE_DIRECTORY = "mapping";

    private static final String BAM_DIFF_EXTENSION = ".prepostrealign.diff";
    private static final String IGNORE_FOR_DIFF_COUNT_PATTERN_2 = ">";
    private static final String IGNORE_FOR_DIFF_COUNT_PATTERN_1 = "<";

    private static final String FLAGSTAT_EXTENSION = ".postrealign.sliced.flagstat";
    private static final String FLAGSTAT_MAPPED_PATTERN = "mapped";
    private static final String FLAGSTAT_END_OF_MAPPED_VALUE_PATTERN = "+";

    public RealignerChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.REALIGNER;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext)
            throws IOException, HealthChecksException {
        final HealthCheck refCheck = extractCheckForSample(runContext.runDirectory(), runContext.refSample());
        final HealthCheck tumorCheck = extractCheckForSample(runContext.runDirectory(), runContext.tumorSample());

        return toPatientResult(refCheck, tumorCheck);
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        return toPatientResult(new HealthCheck(runContext.refSample(),
                                               REALIGNER_CHECK_NAME,
                                               HealthCheckConstants.ERROR_VALUE),
                               new HealthCheck(runContext.tumorSample(),
                                               REALIGNER_CHECK_NAME,
                                               HealthCheckConstants.ERROR_VALUE));
    }

    @NotNull
    private BaseResult toPatientResult(@NotNull final HealthCheck refCheck, @NotNull final HealthCheck tumorCheck) {
        refCheck.log(LOGGER);
        tumorCheck.log(LOGGER);
        return new PatientResult(checkType(),
                                 Collections.singletonList(refCheck),
                                 Collections.singletonList(tumorCheck));
    }

    @NotNull
    private static HealthCheck extractCheckForSample(@NotNull final String runDirectory,
                                                     @NotNull final String sampleId)
            throws IOException, HealthChecksException {
        final String basePath = getBasePathForSample(runDirectory, sampleId);

        final Path bamDiffPath = PathPrefixSuffixFinder.build().findPath(basePath, sampleId, BAM_DIFF_EXTENSION);
        final long diffCount = readDiffCountFromBamDiff(bamDiffPath);

        final Path flagStatPath = PathPrefixSuffixFinder.build().findPath(basePath, sampleId, FLAGSTAT_EXTENSION);
        final long mappedValue = readMappedFromFlagstat(flagStatPath);

        final String value = new DecimalFormat(REALIGNER_CHECK_PRECISION).format((double) diffCount / mappedValue);
        return new HealthCheck(sampleId, REALIGNER_CHECK_NAME, value);
    }

    private static long readMappedFromFlagstat(@NotNull final Path flagStatPath)
            throws IOException, HealthChecksException {
        final List<String> lines = FileReader.build().readLines(flagStatPath);

        final Optional<String> mappedLine = lines.stream()
                                                 .filter(line -> line.contains(FLAGSTAT_MAPPED_PATTERN))
                                                 .findFirst();
        if (!mappedLine.isPresent()) {
            throw new LineNotFoundException(flagStatPath.toString(), FLAGSTAT_MAPPED_PATTERN);
        }
        final String mapped = mappedLine.get();
        if (!mapped.contains(FLAGSTAT_END_OF_MAPPED_VALUE_PATTERN)) {
            throw new MalformedFileException(String.format(MALFORMED_FILE_MSG,
                                                           flagStatPath.toString(),
                                                           FLAGSTAT_END_OF_MAPPED_VALUE_PATTERN));
        }
        final String mappedValue = mapped.substring(0, mapped.indexOf(FLAGSTAT_END_OF_MAPPED_VALUE_PATTERN));
        return Long.valueOf(mappedValue.trim());
    }

    private static long readDiffCountFromBamDiff(@NotNull final Path bamDiffPath)
            throws IOException, HealthChecksException {
        final List<String> lines = FileReader.build().readLines(bamDiffPath);
        return lines.stream()
                    .filter(line -> !line.startsWith(IGNORE_FOR_DIFF_COUNT_PATTERN_1)
                                    && !line.startsWith(IGNORE_FOR_DIFF_COUNT_PATTERN_2))
                    .count();
    }

    @NotNull
    private static String getBasePathForSample(@NotNull final String runDirectory, @NotNull final String sampleId) {
        return runDirectory + File.separator + sampleId + File.separator + REALIGNER_BASE_DIRECTORY;
    }
}
