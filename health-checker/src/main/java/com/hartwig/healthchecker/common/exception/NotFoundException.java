package com.hartwig.healthchecker.common.exception;

import org.jetbrains.annotations.NotNull;

public class NotFoundException extends HealthChecksException {

    private static final long serialVersionUID = 6625914789247983088L;

    public NotFoundException(@NotNull final String message) {
        super(message);
    }
}
