package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PurpleChecker extends ErrorHandlingChecker {

    private static final Logger LOGGER = LogManager.getLogger(PurpleChecker.class);

    @NotNull
    @Override
    protected BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final String purpleDirectory = runContext.runDirectory() + File.separator + "purple";
        final PurpleQC qcCheck = PurpleQCFile.read(PurpleQCFile.generateFilename(purpleDirectory, runContext.tumorSample()));
        final List<HealthCheck> checks = Lists.newArrayList();

        final String segmentScore = String.valueOf(qcCheck.segmentScore());
        checks.add(new HealthCheck(runContext.tumorSample(), PurpleCheck.PURPLE_SEGMENT_SCORE.toString(), segmentScore));
        checks.add(new HealthCheck(runContext.tumorSample(), PurpleCheck.AMBER_GENDER.toString(), qcCheck.amberGender().toString()));
        checks.add(new HealthCheck(runContext.tumorSample(), PurpleCheck.COBALT_GENDER.toString(), qcCheck.cobaltGender().toString()));

        return toMultiValueResult(checks);
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        if (runContext.isSomaticRun()) {
            final List<HealthCheck> checks = Lists.newArrayList();

            checks.add(new HealthCheck(runContext.tumorSample(), PurpleCheck.PURPLE_SEGMENT_SCORE.toString(),
                    HealthCheckConstants.ERROR_VALUE));
            checks.add(new HealthCheck(runContext.tumorSample(), PurpleCheck.AMBER_GENDER.toString(), HealthCheckConstants.ERROR_VALUE));
            checks.add(new HealthCheck(runContext.tumorSample(), PurpleCheck.COBALT_GENDER.toString(), HealthCheckConstants.ERROR_VALUE));

            return toMultiValueResult(checks);
        } else {
            return new NoResult(CheckType.PURPLE);
        }
    }

    @NotNull
    private BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);
        return new MultiValueResult(CheckType.PURPLE, checks);
    }
}
