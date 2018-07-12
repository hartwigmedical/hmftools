package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PurpleChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(PurpleChecker.class);

    @NotNull
    @Override
    public BaseResult run(@NotNull final RunContext runContext) {
        if (!runContext.isSomaticRun()) {
            return new NoResult(CheckType.PURPLE);
        }
        final String purpleDirectory = runContext.runDirectory() + File.separator + "purple";
        String path = PurpleQCFile.generateFilename(purpleDirectory, runContext.tumorSample());
        List<String> lines;
        try {
            lines = LineReader.build().readLines(new File(path).toPath(), x -> x.contains("QCStatus"));
        } catch (IOException exc) {
            LOGGER.warn("Could not load purple qc file.");
            return new NoResult(CheckType.PURPLE);
        }

        assert lines.size() == 1;

        final List<HealthCheck> checks = Lists.newArrayList();
        String[] parts = lines.get(0).split("\t");

        checks.add(new HealthCheck(runContext.tumorSample(), PurpleCheck.QC_STATUS.toString(), parts[1]));

        return toMultiValueResult(checks);
    }

    @NotNull
    private BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);
        return new MultiValueResult(CheckType.PURPLE, checks);
    }
}
