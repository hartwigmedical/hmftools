package com.hartwig.hmftools.common.exception;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public class EmptyFileException extends IOException {

    private static final String EMPTY_FILES_ERROR = "File %s was found empty in path -> %s";

    public EmptyFileException(@NotNull final String fileName, final String filePath) {
        super(String.format(EMPTY_FILES_ERROR, fileName, filePath));
    }
}
