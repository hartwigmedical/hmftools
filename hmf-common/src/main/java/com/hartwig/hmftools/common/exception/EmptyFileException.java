package com.hartwig.hmftools.common.exception;

import org.jetbrains.annotations.NotNull;

public class EmptyFileException extends HartwigException {

    private static final long serialVersionUID = -8396650626359037492L;

    private static final String EMPTY_FILES_ERROR = "File %s was found empty in path -> %s";

    public EmptyFileException(@NotNull final String fileName, final String filePath) {
        super(String.format(EMPTY_FILES_ERROR, fileName, filePath));
    }
}
