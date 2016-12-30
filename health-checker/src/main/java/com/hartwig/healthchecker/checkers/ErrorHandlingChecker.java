package com.hartwig.healthchecker.checkers;

import java.io.IOException;

import com.hartwig.healthchecker.exception.HealthChecksException;
import com.hartwig.healthchecker.io.dir.RunContext;
import com.hartwig.healthchecker.result.BaseResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public abstract class ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(ErrorHandlingChecker.class);
    private static final String EXCEPTION_MSG = "Got an exception with message: %s";

    @NotNull
    @Override
    public BaseResult run(@NotNull final RunContext runContext) {
        BaseResult result;
        try {
            result = tryRun(runContext);
        } catch (IOException | HealthChecksException exception) {
            LOGGER.error(String.format(EXCEPTION_MSG, exception.getMessage()));
            result = errorRun(runContext);
        }
        return result;
    }

    @NotNull
    protected abstract BaseResult tryRun(@NotNull RunContext runContext)
            throws IOException, HealthChecksException;

    @NotNull
    protected abstract BaseResult errorRun(@NotNull RunContext runContext);
}
