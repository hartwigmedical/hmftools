package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;
import com.hartwig.hmftools.healthchecker.result.SingleValueResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class AmberChecker extends ErrorHandlingChecker {

    private static final Logger LOGGER = LogManager.getLogger(AmberChecker.class);

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.AMBER;
    }

    @NotNull
    @Override
    protected BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final String amberDirectory = runContext.runDirectory() + File.separator + "amber";
        final AmberQC qcCheck = AmberQCFile.read(AmberQCFile.generateFilename(amberDirectory, runContext.tumorSample()));

        final String meanBaf = String.valueOf(qcCheck.meanBAF());
        final HealthCheck healthCheck = new HealthCheck(runContext.tumorSample(), AmberCheck.MEAN_BAF.toString(), meanBaf);
        return toSingleValueResult(healthCheck);
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        if (runContext.isSomaticRun()) {
            final HealthCheck healthCheck =
                    new HealthCheck(runContext.tumorSample(), AmberCheck.MEAN_BAF.toString(), HealthCheckConstants.ERROR_VALUE);
            return toSingleValueResult(healthCheck);
        } else {
            return new NoResult(checkType());
        }
    }

    @NotNull
    private BaseResult toSingleValueResult(@NotNull final HealthCheck check) {
        check.log(LOGGER);
        return new SingleValueResult(checkType(), check);
    }
}
