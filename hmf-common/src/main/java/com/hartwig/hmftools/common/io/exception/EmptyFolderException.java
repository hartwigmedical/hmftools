package com.hartwig.hmftools.common.io.exception;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public class EmptyFolderException extends IOException {

    private static final String MESSAGE = "Folder %s is empty";

    public EmptyFolderException(@NotNull final String folder) {
        super(String.format(MESSAGE, folder));
    }
}
