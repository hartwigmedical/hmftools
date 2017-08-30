package com.hartwig.hmftools.healthchecker.runners;

import java.util.List;

import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HealthCheck {

    private static final String LOG_MSG = "Check '%s' for sample '%s' has value '%s'";

    @NotNull
    private final String sampleId;
    @NotNull
    private final String checkName;
    @NotNull
    private final String value;

    HealthCheck(@NotNull final String sampleId, @NotNull final String checkName,
            @NotNull final String value) {
        this.sampleId = sampleId;
        this.checkName = checkName;
        this.value = value;
    }

    @NotNull
    String getSampleId() {
        return sampleId;
    }

    @NotNull
    String getCheckName() {
        return checkName;
    }

    @NotNull
    String getValue() {
        return value;
    }

    void log(@NotNull Logger logger) {
        logger.info(String.format(LOG_MSG, checkName, sampleId, value));
    }

    static void log(@NotNull final Logger logger, @NotNull final List<HealthCheck> reports) {
        for (final HealthCheck report : reports) {
            report.log(logger);
        }
    }
}
