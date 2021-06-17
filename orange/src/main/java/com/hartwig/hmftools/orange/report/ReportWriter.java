package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(ReportWriter.class);

    @NotNull
    private final String outputDir;

    public ReportWriter(@NotNull final String outputDir) {
        this.outputDir = outputDir;
    }

    public void write(@NotNull OrangeReport report) {
        LOGGER.info("Writing {}", report);
    }
}
