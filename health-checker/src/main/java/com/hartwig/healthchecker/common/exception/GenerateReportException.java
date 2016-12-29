package com.hartwig.healthchecker.common.exception;

import org.jetbrains.annotations.NotNull;

public class GenerateReportException extends HealthChecksException {

    private static final long serialVersionUID = 7343366902433673704L;

    public GenerateReportException(@NotNull final String message) {
        super(message);
    }
}
