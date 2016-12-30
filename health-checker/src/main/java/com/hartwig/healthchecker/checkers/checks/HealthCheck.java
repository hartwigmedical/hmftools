package com.hartwig.healthchecker.checkers.checks;

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

    public HealthCheck(@NotNull final String sampleId, @NotNull final String checkName,
            @NotNull final String value) {
        this.sampleId = sampleId;
        this.checkName = checkName;
        this.value = value;
    }

    @NotNull
    public String getSampleId() {
        return sampleId;
    }

    @NotNull
    public String getCheckName() {
        return checkName;
    }

    @NotNull
    public String getValue() {
        return value;
    }

    public void log(@NotNull Logger logger) {
        logger.info(String.format(LOG_MSG, checkName, sampleId, value));
    }

    public static void log(@NotNull Logger logger, @NotNull List<HealthCheck> reports) {
        for (HealthCheck report : reports) {
            report.log(logger);
        }
    }
}
