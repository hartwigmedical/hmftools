package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.exception.MalformedFileException;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.SingleValueResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.KinshipCheck;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.KINSHIP)
public class KinshipChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(KinshipChecker.class);

    private static final String MALFORMED_FILE_MSG = "Malformed %s file is path %s -> %s lines found was expecting %s";

    private static final String KINSHIP_EXTENSION = ".kinship";
    private static final int EXPECTED_NUM_LINES = 2;
    private static final String COLUMN_SEPARATOR = "\t";
    private static final int KINSHIP_COLUMN = 7;

    public KinshipChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.KINSHIP;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final Path kinshipPath = PathExtensionFinder.build().findPath(runContext.runDirectory(), KINSHIP_EXTENSION);
        final List<String> kinshipLines = FileReader.build().readLines(kinshipPath);
        if (kinshipLines.size() != EXPECTED_NUM_LINES) {
            throw new MalformedFileException(
                    String.format(MALFORMED_FILE_MSG, KinshipCheck.KINSHIP_TEST.toString(), runContext.runDirectory(),
                            kinshipLines.size(), EXPECTED_NUM_LINES));
        }
        final Optional<HealthCheck> optCheck = kinshipLines.stream().skip(1).map(line -> {
            final String[] values = line.split(COLUMN_SEPARATOR);
            return new HealthCheck(runContext.tumorSample(), KinshipCheck.KINSHIP_TEST.toString(),
                    values[KINSHIP_COLUMN]);
        }).findFirst();

        assert optCheck.isPresent();

        return toSingleValueResult(optCheck.get());
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        return toSingleValueResult(new HealthCheck(runContext.tumorSample(), KinshipCheck.KINSHIP_TEST.toString(),
                HealthCheckConstants.ERROR_VALUE));
    }

    @NotNull
    private BaseResult toSingleValueResult(@NotNull final HealthCheck check) {
        check.log(LOGGER);
        return new SingleValueResult(checkType(), check);
    }
}
