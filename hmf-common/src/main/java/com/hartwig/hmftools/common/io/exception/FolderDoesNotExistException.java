package com.hartwig.hmftools.common.io.exception;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public class FolderDoesNotExistException extends IOException {

    private static final String MESSAGE = "Folder %s does not exist";

    public FolderDoesNotExistException(@NotNull final String folderName) {
        super(String.format(MESSAGE, folderName));
    }
}
