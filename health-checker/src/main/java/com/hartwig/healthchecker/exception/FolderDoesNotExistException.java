package com.hartwig.healthchecker.exception;

import org.jetbrains.annotations.NotNull;

public class FolderDoesNotExistException extends HealthChecksException {

    private static final long serialVersionUID = -8396650626359037492L;

    private static final String MESSAGE = "Folder %s does not exist";

    public FolderDoesNotExistException(@NotNull final String folderName) {
        super(String.format(MESSAGE, folderName));
    }
}
