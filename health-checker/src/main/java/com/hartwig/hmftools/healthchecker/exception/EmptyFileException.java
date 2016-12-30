package com.hartwig.hmftools.healthchecker.exception;

import org.jetbrains.annotations.NotNull;

public class EmptyFileException extends HealthChecksException {

    private static final long serialVersionUID = -8396650626359037492L;

    private static final String EMPTY_FILES_ERROR = "File %s was found empty in path -> %s";

    public EmptyFileException(@NotNull final String fileName, final String filePath) {
        super(String.format(EMPTY_FILES_ERROR, fileName, filePath));
    }
}
