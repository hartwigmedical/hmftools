package com.hartwig.healthchecker.common.exception;

import org.jetbrains.annotations.NotNull;

public class HealthChecksException extends Exception {

    private static final long serialVersionUID = -7416221606300503967L;

    public HealthChecksException(@NotNull final String message) {
        super(message);
    }
}
