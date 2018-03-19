package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;
import com.hartwig.hmftools.healthchecker.result.SingleValueResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class AmberChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(AmberChecker.class);

    @NotNull
    public BaseResult run(@NotNull final RunContext runContext) throws IOException {
        if (!runContext.isSomaticRun()) {
            return new NoResult(CheckType.AMBER);
        }

        final String amberDirectory = runContext.runDirectory() + File.separator + "amber";
        final AmberQC qcCheck = AmberQCFile.read(AmberQCFile.generateFilename(amberDirectory, runContext.tumorSample()));

        final String meanBaf = String.valueOf(qcCheck.meanBAF());
        final HealthCheck healthCheck = new HealthCheck(runContext.tumorSample(), AmberCheck.MEAN_BAF.toString(), meanBaf);
        return toSingleValueResult(healthCheck);
    }

    @NotNull
    private static BaseResult toSingleValueResult(@NotNull final HealthCheck check) {
        check.log(LOGGER);
        return new SingleValueResult(CheckType.AMBER, check);
    }
}
