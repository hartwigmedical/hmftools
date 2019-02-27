package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class AmberChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(AmberChecker.class);

    @NotNull
    @Override
    public BaseResult run(@NotNull final RunContext runContext) {
        if (!runContext.isSomaticRun()) {
            return new NoResult(CheckType.AMBER);
        }

        final AmberQC qcCheck;
        try {
            final String amberDirectory = runContext.runDirectory() + File.separator + "amber";
            qcCheck = AmberQCFile.read(AmberQCFile.generateFilename(amberDirectory, runContext.tumorSample()));
        } catch (IOException exception) {
            LOGGER.warn("Could not load amber qc file.");
            return new NoResult(CheckType.AMBER);
        }

        final HealthCheck bafCheck = new HealthCheck(runContext.tumorSample(), AmberCheck.MEAN_BAF.toString(), String.valueOf(qcCheck.meanBAF()));
        final HealthCheck contaminationCheck = new HealthCheck(runContext.tumorSample(), AmberCheck.CONTAMINATION.toString(), String.valueOf(qcCheck.contamination()));
        return toMultiValueResult(Lists.newArrayList(bafCheck, contaminationCheck));
    }

    @NotNull
    private static BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);

        return new MultiValueResult(CheckType.AMBER, checks);
    }
}
